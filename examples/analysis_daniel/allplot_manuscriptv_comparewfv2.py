"""Created on Wed Mar  6 12:35:05 2024.

@author: dpm42
"""

import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import tol_colors as tc
from scipy import stats
from scipy.stats import pearsonr

os.chdir("../../")
import sys

sys.path.pop(-2)
from src.core.plotter import datamatch_merge

mpl.rcParams["figure.dpi"] = 400
pd.options.mode.chained_assignment = None


def concordance_correlation_coefficient(y_true, y_pred):
    """Concordance correlation coefficient."""
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    # Raw data
    dct = {"y_true": y_true, "y_pred": y_pred}
    df = pd.DataFrame(dct)
    # Remove NaNs
    df = df.dropna()
    # Pearson product-moment correlation coefficients
    y_true = df["y_true"]
    y_pred = df["y_pred"]
    cor = np.corrcoef(y_true, y_pred)[0][1]
    # Means
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)
    # Population variances
    var_true = np.var(y_true)
    var_pred = np.var(y_pred)
    # Population standard deviations
    sd_true = np.std(y_true)
    sd_pred = np.std(y_pred)
    # Calculate CCC
    numerator = 2 * cor * sd_true * sd_pred
    denominator = var_true + var_pred + (mean_true - mean_pred) ** 2

    return numerator / denominator


def CCC_wp(y_true, y_pred):
    """Concordance correlation coefficient.

    With pval
    """
    r = concordance_correlation_coefficient(y_true, y_pred)
    t = r * np.sqrt(len(y_true) - 2) / np.sqrt(1 - r**2)
    p = 2 * (1 - stats.t.cdf(abs(t), len(y_true) - 2))
    return r, p


def corrcalc(data, comparison):
    corrs = []
    for nsim in pd.unique(data["nsim"]):
        # ax.set_title(f'fiber diam: {s}μm')
        corr = {}
        for sample in pd.unique(data["sample"]):
            thisdat = data[(data["nsim"] == nsim) & (data["sample"] == sample)]
            corr[sample] = round(pearsonr(thisdat[comparison[0]], thisdat[comparison[1]])[0], 3)
        corrs.append(corr)
    return corrs


def addpwfd(data, sim, infile="plotconfig"):
    with open(f"examples/analysis_daniel/{infile}.json") as f:
        config = json.load(f)
    nsim_key = config["sim_data"][sim]["nsim_key"]
    for nsim in nsim_key:
        pulse_width = nsim_key[nsim]["pulse_width"]
        fiber_diam = nsim_key[nsim]["fiber_diam"]
        nsim = int(nsim)
        data.loc[data.nsim == nsim, "pulse_width"] = pulse_width
        data.loc[data.nsim == nsim, "fiber_diam"] = fiber_diam
    return data


def pe(correct, est):
    """Calculate the percent error.

    :param correct: correct value
    :param est: estimated value
    :return: percent error
    """
    return 100 * abs(est - correct) / correct


def pe_noabs(correct, est, doabs=False):
    """Calculate the percent error.

    :param correct: correct value
    :param est: estimated value
    :return: percent error
    """
    diff = est - correct
    if doabs:
        diff = np.abs(diff)
    return 100 * diff / correct


import math


def compute_reorder_cost(order1, order2):
    assert set(order1) == set(order2)
    assert len(set(order1)) == len(order1)
    assert len(set(order2)) == len(order2)
    sumdist = 0
    for item in order1:
        distance = abs(order1.index(item) - order2.index(item))
        sumdist += distance
    if len(order1) % 2 == 0:
        maxlen = sum(range(len(order1))[int(math.ceil(0.5 * len(order1))) :]) * 2 / len(order1)
    else:
        maxlen = (
            range(len(order1))[int(math.floor(0.5 * len(order1)))]
            + (sum(range(len(order1))[int(math.ceil(0.5 * len(order1))) :])) * 2
        ) / len(order1)
    return sumdist / len(order1) / maxlen


def calculate_dose_response(df, threshold_column, outcolumn, grouping_columns):
    def calculate_dose_response_per_sample(group):
        threshold_values = group[threshold_column].values
        activated_fibers = group[threshold_column].values[:, None] <= threshold_values
        dose_response = activated_fibers.sum(axis=0) / len(group)
        group[outcolumn] = dose_response
        return group

    return df.groupby(grouping_columns).apply(calculate_dose_response_per_sample)


def rownames(g, row_template, **ylabel_kws):
    args = dict(row_var=g._row_var)
    # Draw the row titles on the left edge of the grid
    for i, row_name in enumerate(g.row_names):
        ax = g.axes[i, 0]
        args.update(dict(row_name=row_name))
        title = row_template.format(**args)
        ax.set_ylabel(title, **ylabel_kws)


# %%
newdefdat_allwf = []
deftomatch = []
allconcats = []
for simNUM, stimtype in zip(
    ['333', '330', '3', '3'],
    ['Asymmetric & Bipolar', 'Asymmetric & Monopolar', 'Symmetric & Bipolar', 'Symmetric & Bipolar (Select)'],
):

    gogo = "initial"

    cath_comparison = ["cathodic", "3D", 3]
    center_comparison = ["center", "3D", 3]
    an_comparison = ["anodic", "3D", 3]
    comparisons = [cath_comparison, center_comparison, an_comparison]
    main_comparison = cath_comparison  # NOTE: changed from center

    pal2d3d = ["#d95f02", "#7570b3"]

    # code whole file to optionally run on 3 or 10 so that I can run all monpolar data
    # %% BEGIN DEFORMATION ANALYSIS
    print("deformation")
    # match fascicle thresholds
    sns.set(font_scale=1, style='white', context='paper')
    threshes = pd.concat(
        [
            pd.read_csv(f"thresh_unmatched_sim{simNUM}_og.csv"),
            pd.read_csv(f"thresh_unmatched_sim{simNUM}_def.csv"),
        ]
    )
    threshes["waveform"] = stimtype
    # remove all rows where nerve label contains "asc" and "sample" contains "3"
    threshes = threshes[
        ~((threshes["nerve_label"].str.contains("asc")) & (threshes["sample"].astype(str).str[2] == "3"))
    ]
    # duplicate all samples where nerve label contains "def" and "sample" contains "3"
    defdat = threshes[(threshes["nerve_label"].str.contains("def")) & (threshes["sample"].astype(str).str[2] == "3")]
    # replacee all "1" with "9" in sample column
    defdat["sample"] = defdat["sample"].astype(str).str.replace("1", "9").astype(int)
    # replace all "def" with "asc" in nerve label column
    defdat["nerve_label"] = defdat["nerve_label"].str.replace("def", "asc")
    # combine with original data
    threshes = pd.concat([threshes, defdat])
    threshes["deformation"] = None
    # where nerve label contains "def" deformation is "Structural"
    threshes.loc[threshes["nerve_label"].str.contains("def"), "deformation"] = "Structural"
    # where nerve label contains "asc" deformation is "ASCENT"
    threshes.loc[threshes["nerve_label"].str.contains("asc"), "deformation"] = "ASCENT"
    # else deformation is "none"
    threshes.loc[threshes["deformation"].isna(), "deformation"] = "Undeformed"
    # strip "def" and "asc" from nerve labels
    threshes["nerve_label"] = threshes["nerve_label"].str.replace("def", "").str.replace("asc", "")

    newdefdat = threshes.copy()

    # remove all nerve_label = '6L'
    newdefdat = newdefdat[newdefdat["nerve_label"] != "6L"]
    newdefdat = newdefdat[newdefdat["nerve_label"] != "2R"]

    newdefdat["type"] = newdefdat["type"].replace({"2D": "extrusion", "3D": "true-3D"})
    newdefdat = addpwfd(newdefdat, simNUM, infile="plotconfig_og")
    newdefdat["fiber_diam"] = newdefdat["fiber_diam"].astype(int)
    # contact is cathodic if the third digit of sample int is 2, anodic if 0
    newdefdat["contact"] = (
        newdefdat["sample"].astype(str).str[2].replace({"2": "cathodic", "1": "center", "0": "anodic", "3": "3D"})
    )

    # set deformation as ordered categorical
    newdefdat["deformation"] = pd.Categorical(
        newdefdat["deformation"],
        categories=["Undeformed", "ASCENT", "Structural"],
        ordered=True,
    )
    newdefdat["nerve_label"] = pd.Categorical(
        newdefdat["nerve_label"], categories=["2L", "3R", "5R", "6R"], ordered=True
    )
    newdefdat["deformed"] = newdefdat["deformation"] != "Undeformed"
    newdefdat.reset_index(inplace=True)

    if 'Select' in stimtype:
        # Filter out cathodic and anodic data
        cathodic_df = newdefdat[newdefdat['contact'] == 'cathodic']
        anodic_df = newdefdat[newdefdat['contact'] == 'anodic']

        # Perform a merge on the relevant columns to get corresponding anodic thresholds for cathodic rows
        merged_df = cathodic_df.merge(
            anodic_df,
            on=['nerve_label', 'nsim', 'waveform', 'deformation', 'master_fiber_index'],
            suffixes=('_cathodic', '_anodic'),
        )

        # Identify rows where cathodic threshold is greater than anodic threshold
        condition = merged_df['threshold_cathodic'] > merged_df['threshold_anodic']

        # Update the thresholds in the original DataFrame using the indices of the cathodic rows
        cathodic_indices_to_update = merged_df.loc[condition, 'threshold_cathodic'].index

        for idx in cathodic_indices_to_update:
            newdefdat.loc[cathodic_df.index[idx], 'threshold'] = merged_df.loc[idx, 'threshold_anodic']

    # all wf data
    newdefdat_allwf.append(newdefdat)

    deftomatch_this = newdefdat.copy()
    deftomatch.append(deftomatch_this)

    # remove unused colors from palette
    defpal = [sns.color_palette("colorblind")[ind] for ind in [0, 2, 3, 5]]
    defdefcomp = newdefdat.query('type=="true-3D" and deformation != "ASCENT"')
    defdefcomp["deformed"] = defdefcomp["deformation"] != "Undeformed"
    defsamples = ["2L", "3R", "5R", "6R"]
    # %% inners need to match the cathodic leading contact
    definnermatch = deftomatch_this.query('deformation!="ASCENT"').reset_index(drop=True)

    # Query the rows with type '3D'
    df_3d = definnermatch.query("type == 'true-3D'")

    # Query the rows with type '2D' and sample ending in '1'
    df_2d = definnermatch.query(f"type == 'extrusion' and contact in {main_comparison}")

    # Merge the 3D and 2D data, keeping track of original row indices
    merged_df = pd.merge(
        df_3d,
        df_2d,
        on=["nerve_label", "master_fiber_index", "nsim", "deformation"],
        suffixes=("_3d", "_2d"),
        how="left",
    )  # TODO remove this how

    # Update the 'inner', 'outer', and 'fiber' columns in the original DataFrame
    definnermatch.loc[df_3d.index, "inner"] = merged_df["inner_2d"].values
    definnermatch.loc[df_3d.index, "outer"] = merged_df["outer_2d"].values
    definnermatch.loc[df_3d.index, "fiber"] = merged_df["fiber_2d"].values
    # %% new idea analysis
    concats = []
    for deftype in ["Undeformed", "Structural"]:
        thismatch = deftomatch_this.query("deformation==@deftype")
        matchednow = datamatch_merge(
            thismatch.query('type=="extrusion"'),
            thismatch.query('type=="true-3D"'),
            "threshold",
            merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index", "waveform"],
        ).drop(columns="type")
        matchednow["deformation"] = deftype
        concats.append(matchednow)
    concats = pd.concat(concats)
    concats["deformation"] = pd.Categorical(
        concats["deformation"], categories=["Undeformed", "Structural"], ordered=True
    )
    concats['waveform'] = stimtype
    allconcats.append(concats)

    # %% dose-response
    print("deformdr")
    defdr = newdefdat.copy().rename(columns={"sample": "samplenum", "type": "modeltype"})
    defdr = calculate_dose_response(
        defdr,
        "threshold",
        "percent_activated",
        grouping_columns=[
            "modeltype",
            "samplenum",
            "fiber_diam",
            "sim",
            "deformation",
            "pulse_width",
        ],
    )
    defdrref = calculate_dose_response(
        defdr,
        "threshold",
        "percent_activated",
        grouping_columns=[
            "modeltype",
            "samplenum",
            "fiber_diam",
            "sim",
            "deformation",
            "pulse_width",
        ],
    )
    import matplotlib.colors as colors

    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
            cmap(np.linspace(minval, maxval, n)),
        )
        return new_cmap

    rdup = truncate_colormap(plt.get_cmap("RdPu"), minval=0.25)
    # %%Setup def
    print("matcheddef")
    # true-3D and extrusion trhresholds as different col same row
    matched_deformation = datamatch_merge(
        deftomatch_this.query('type=="extrusion" and deformation!="ASCENT"'),
        deftomatch_this.query('type=="true-3D" and deformation!="ASCENT"'),
        "threshold",
        merge_cols=[
            "model",
            "sim",
            "nerve_label",
            "nsim",
            "master_fiber_index",
            "deformation",
        ],
    ).drop(columns="type")
    # for matched boxplots
    repeated_deformation = matched_deformation.melt(
        id_vars=[x for x in matched_deformation.columns if "threshold" not in x],
        value_vars=["threshold3d", "threshold"],
        var_name="type",
        value_name="threshold",
    )
    repeated_deformation["type"] = repeated_deformation["type"].replace(
        {"threshold": "extrusion", "threshold3d": "true-3D"}
    )
    repeated_deformation["type"] = pd.Categorical(
        repeated_deformation["type"], categories=["extrusion", "true-3D"], ordered=True
    )
# %%
newdefdat_allwf = pd.concat(newdefdat_allwf).query('fiber_diam in [3,13]')
deftomatch = pd.concat(deftomatch)
allconcats = pd.concat(allconcats)
sns.set(context='paper', style='white')
sys.exit("prepdone")
# %% compare dose response across waveforms and deformation
peses = []
pemeans = []
onsets_sats = {}
for comparison in comparisons:
    for stringdat in ["Undeformed", "Structural", "ASCENT"]:
        for waveform in newdefdat_allwf.waveform.unique():
            thiscontact = comparison[0]
            subdat = newdefdat_allwf.query(
                f"deformation=='{stringdat}' and contact in {comparison} and waveform=='{waveform}'"
            )
            levels = {
                "onset": 10,
                "half": 50,
                "saturation": 90,
            }
            grouped = subdat.groupby(
                [
                    "sample",
                    "fiber_diam",
                    "type",
                    "sim",
                    "nerve_label",
                    "model",
                    "nsim",
                    "deformation",
                ]
            )
            analysis = grouped.agg(
                {
                    "threshold": [
                        lambda x: np.percentile(x, q=levels["onset"]),
                        lambda x: np.percentile(x, q=levels["half"]),
                        lambda x: np.percentile(x, q=levels["saturation"]),
                    ]
                }
            )
            analysis.columns = ["_".join(col_name).rstrip("_") for col_name in analysis.columns]
            analysis.rename(
                columns={
                    "threshold_<lambda_0>": "onset",
                    "threshold_<lambda_1>": "half",
                    "threshold_<lambda_2>": "saturation",
                },
                inplace=True,
            )
            analysis = analysis.reset_index()
            # combine onset, saturation, and half into one column with identifier
            compiled_data = analysis.melt(
                id_vars=[
                    "sample",
                    "fiber_diam",
                    "sim",
                    "type",
                    "nerve_label",
                    "model",
                    "nsim",
                ],
                value_vars=["onset", "half", "saturation"],
                var_name="level",
                value_name="threshold",
            )

            # set up facetgrid with nsim as row and level as columns
            compiled_data.reset_index(inplace=True)
            # set fiber_diam to category
            compiled_data.type = compiled_data.type.astype("category")
            # add a units column with unique number for each combination of fiber_diam and level
            compiled_data["units"] = compiled_data.groupby(["fiber_diam", "level", "nerve_label"]).ngroup()
            compiled_data["fiber_diam"] = compiled_data["fiber_diam"].astype(int)
            compiled_data.dropna(inplace=True)

            # calculate percent error for sample onset and saturation as well as population onset and saturation
            pes = []
            for level in ["onset", "half", "saturation"]:
                for sam in compiled_data.nerve_label.unique():
                    for fiber_diam in compiled_data.fiber_diam.unique():
                        a = compiled_data.query(
                            f"nerve_label == '{sam}' and level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam}"
                        )["threshold"].values
                        b = compiled_data.query(
                            f"nerve_label == '{sam}' and level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}"
                        )["threshold"].values
                        assert len(a) == len(b) == 1
                        pe_res = pe_noabs(b[0], a[0], doabs=False)
                        pes.append(
                            {
                                "level": level,
                                "nerve_label": sam,
                                "fiber_diam": fiber_diam,
                                "pe": pe_res,
                            }
                        )
            pes = pd.DataFrame(pes)
            pes["deformation"] = stringdat
            pes["contact"] = thiscontact
            pes["wf"] = waveform
            assert thiscontact != np.nan
            peses.append(pes)
            print("Max 3 um", stringdat, np.amax(pes.query("fiber_diam==3").pe))
            print("Max 13 um", stringdat, np.amax(pes.query("fiber_diam==13").pe))

            # now calculate percent error for population onset and saturation
            pemean = []
            for level in ["onset", "half", "saturation"]:
                for fiber_diam in compiled_data.fiber_diam.unique():
                    a = compiled_data.query(
                        f"level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam}"
                    )["threshold"].values
                    b = compiled_data.query(f"level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}")[
                        "threshold"
                    ].values
                    pe_res = pe_noabs(np.median(b), np.median(a), doabs=False)
                    pemean.append({"level": level, "fiber_diam": fiber_diam, "pe": pe_res})

            pemean = pd.DataFrame(pemean)
            pemean["deformation"] = stringdat
            pemean["contact"] = thiscontact
            pemean["wf"] = waveform
            pemeans.append(pemean)
            print("Max 3 um median", stringdat, np.amax(pemean.query("fiber_diam==3").pe))
            print("Max 13 um median", stringdat, np.amax(pemean.query("fiber_diam==13").pe))

            onsets_sats[stringdat] = compiled_data

sns.set(style="white", context="paper")
allpes = pd.concat(peses)
# allpes['level'].replace({
#     'onset':'10','half':'50','saturation':'90'},inplace=True)
# allpes['deformation'].replace({
#     'ASCENT':'ASCENT','Structural':'Structural','Undeformed':'Undef.'},inplace=True)
# allpes['comb'] = allpes['level'].astype(str) + '\n' + allpes['deformation']
# allpes.sort_values(by=['deformation','level'])
allpes["level"] = pd.Categorical(allpes.level, categories=["onset", "half", "saturation"], ordered=True)
allpes = allpes.query('contact=="cathodic" and deformation=="Structural"')

allpes.sort_values(by=["contact", "deformation"], inplace=True)
allpes["comb"] = allpes["deformation"].astype(str) + "\n" + allpes["contact"].astype(str)

sns.set(context='paper', style='white', font_scale=1)
allpes = allpes.query('deformation!="ASCENT"')
allpes["deformation"] = allpes.deformation.replace({'Structural': 'Deformed'})

allpes.wf = allpes.wf.str.replace(' ', '\n')
allpes.fiber_diam = allpes.fiber_diam.astype(str)
wfpal = [tc.tol_cset('muted')[i] for i in [0, 2, 1, 3]]
# %%
sns.set(context='paper', style='white', font_scale=1)
allpes['pe_abs'] = allpes['pe'].abs()
# rewrite the above with map dataframe
g = sns.FacetGrid(
    data=allpes,
    col="fiber_diam",
    # row_order=["saturation", "half", "onset"],
)
g.map_dataframe(
    sns.pointplot,
    x="level",
    y="pe_abs",
    marker="s",
    estimator="median",
    errorbar=None,
    dodge=0.5,
    hue='wf',
    # hue_order = ['Asymmetric\n&\nBipolar', 'Asymmetric\n&\nMonopolar', 'Symmetric\n&\nBipolar', 'Symmetric\n&\nBipolar\n(Select)'],
    palette=wfpal,
)
g.map_dataframe(
    sns.swarmplot,
    x="level",
    y="pe_abs",
    # marker="s",
    # estimator="median",
    # errorbar=None,
    dodge=0.2,
    hue='wf',
    # hue_order = ['Asymmetric\n&\nBipolar', 'Asymmetric\n&\nMonopolar', 'Symmetric\n&\nBipolar', 'Symmetric\n&\nBipolar\n(Select)'],
    palette=wfpal,
    edgecolor='lightgray',
    linewidth=0.5,
    s=4,
)
plt.gcf().set_size_inches(4, 2.5)
g.set_titles(col_template="{col_name}", row_template="")
rownames(g, row_template="Absolute Percent Difference (%)\n{row_name}")
for ax in g.axes.ravel():
    plt.sca(ax)
    # plt.xticks(rotation=20)
    plt.xlabel("")
    ax.set_xticklabels(['onset', 'half', 'saturation'])
plt.subplots_adjust(hspace=0.1, wspace=0)
plt.gcf().set_size_inches(4, 2)
# plt.ylim(None, 80)
plt.xlim(-0.5, 2.5)
g.add_legend(
    title='',
    ncol=2,
    # labels=['Asymmetric & Bipolar', 'Asymmetric & Monopolar', 'Symmetric & Bipolar', 'Symmetric & Bipolar (Select)'],
)  # ,ncol=2,bbox_to_anchor=[0.6,0.7], frameon=True,framealpha=1)
g.set_ylabels('Absolute Percent Difference (%)')
g.set_titles(col_template='D: {col_name} μm')

# %% recruitment cost:
sns.set(font_scale=1, context='paper', style='white')
sns.set_style("white")
datahere = deftomatch.query('deformation in ["Structural","Undeformed"]')
mathere = allconcats.copy()
scores = []
FASRS = []
for comp in [center_comparison, cath_comparison, an_comparison]:
    for waveform in datahere.waveform.unique():
        for deformation in datahere.deformation.unique():
            threshdat = datahere.query(f"contact in {comp} and deformation==@deformation and waveform==@waveform")
            for nerve in pd.unique(threshdat["nerve_label"]):
                for n in [3, 13]:
                    shortdat = threshdat.query(f'nerve_label=="{nerve}" and fiber_diam=={n}')
                    shortdat['percent_activated'] = np.nan
                    assert len(shortdat) > 0
                    data2d = shortdat.query('type=="extrusion"').sort_values("threshold").master_fiber_index
                    data3d = shortdat.query('type=="true-3D"').sort_values("threshold").master_fiber_index
                    rc = compute_reorder_cost(list(data2d), list(data3d))
                    concathere = mathere.query(
                        f'nerve_label=="{nerve}" and fiber_diam=={n} and contact in {comp} and deformation==@deformation and waveform==@waveform'
                    )
                    assert len(concathere) > 0

                    data2d = concathere.threshold
                    data3d = concathere.threshold3d
                    ccc = concordance_correlation_coefficient(list(data2d), list(data3d))
                    data2d = calculate_dose_response(
                        shortdat.query('type=="extrusion"').sort_values("master_fiber_index"),
                        'threshold',
                        'percent_activated',
                        grouping_columns=['deformation'],
                    ).percent_activated

                    data3d = calculate_dose_response(
                        shortdat.query('type=="true-3D"').sort_values("master_fiber_index"),
                        'threshold',
                        'percent_activated',
                        grouping_columns=['deformation'],
                    ).percent_activated
                    accc = concordance_correlation_coefficient(list(data2d), list(data3d))
                    scores.append(
                        {
                            "sample": nerve,
                            "fiber_diam": n,
                            "score2d3d": rc,
                            "deformation": deformation,
                            "slice": comp[0],
                            "CCC": ccc,
                            "wf": waveform,
                            "aCCC": accc,
                        }
                    )
scoredat = pd.DataFrame(scores)
scoredat["slice"] = scoredat["slice"].replace({"cathodic": "cathodic", "anodic": "anodic", "center": "center"})
scoredat["slice"] = pd.Categorical(scoredat["slice"], categories=["cathodic", "center", "anodic"], ordered=True)
scoredat["fiber_diam"] = pd.Categorical(scoredat["fiber_diam"].astype(int), categories=[3, 13], ordered=True)
scoredat = scoredat.query('slice=="cathodic" and deformation == "Structural"')
# %%
# barplot
g = sns.FacetGrid(data=scoredat, margin_titles=True)
g.map_dataframe(
    sns.barplot,
    y="aCCC",
    x="fiber_diam",
    hue="wf",
    palette=wfpal,
    errorbar=None,
    estimator='median',
    legend=False,
)
g.map_dataframe(
    sns.stripplot,
    y="aCCC",
    x="fiber_diam",
    hue="wf",
    dodge=True,
    color="black",
    edgecolor="white",
    linewidth=0.5,
    legend=False,
)
# plt.axhline(0.45, color='red', ls='--')
# plt.text(0, 0.47, 'random', color='red')

# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel("Activation Reordering")
g.set_titles(col_template="D: {col_name} μm")
g.set_ylabels("aCCC\n(extrusion vs true-3D)")
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
plt.ylim(0, 1)
print(scoredat.groupby(["deformation"]).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
# g.fig.set_size_inches(8,4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel("Fiber Diameter (μm)")
handles, labs = plt.gca().get_legend_handles_labels()

plt.gcf().set_size_inches(2, 2)
g.add_legend(title="", handles=handles[:3], labels=labs[:3])
# second barplot
# barplot
g = sns.FacetGrid(data=scoredat, margin_titles=True)
g.map_dataframe(
    sns.barplot,
    y="CCC",
    x="fiber_diam",
    hue="wf",
    palette=wfpal,
    errorbar=None,
    estimator='median',
    legend=False,
)
g.map_dataframe(
    sns.stripplot,
    y="CCC",
    x="fiber_diam",
    hue="wf",
    dodge=True,
    color="black",
    edgecolor="white",
    linewidth=0.5,
    legend=False,
)
plt.ylabel('tCCC')
# ax.set_xlabel('Fiber Diameter (μm)')
# plt.ylabel("Activation Reordering")
g.set_titles(col_template="D: {col_name} μm")
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
print(scoredat.groupby(["deformation"]).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel("Fiber Diameter (μm)")
handles, labs = plt.gca().get_legend_handles_labels()
g.add_legend(title="", handles=handles[:3], labels=labs[:3])
plt.gcf().set_size_inches(2, 2)

# %% activation position
for waveform in newdefdat_allwf.waveform.unique():
    print(waveform)
    plt.figure()
    sns.set(font_scale=1, style="white", context='paper')
    newthreshz = newdefdat_allwf.query("deformation=='Structural'")
    newthreshz["activation_zpos"] = newthreshz["activation_zpos"] / 10000
    fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
    for nerve_label in pd.unique(newthreshz.nerve_label):
        for ax, modeltype in zip(axs, ["extrusion", "true-3D"]):
            g = sns.histplot(
                data=newthreshz.query(
                    f"fiber_diam in [3,13] and waveform==@waveform and contact in {main_comparison} and nerve_label=='{nerve_label}' and type=='{modeltype}'"
                ).rename(columns={"nerve_label": "Sample"}),
                y="activation_zpos",
                hue="fiber_diam",
                # hue='Sample',
                # facet_kws={'sharex': False},
                # kind='kde',
                palette=[sns.color_palette("binary")[5], sns.color_palette("binary")[2]],
                common_norm=False,
                # legend=False,
                # multiple="fill",
                element="poly",
                fill=False,
                bins=np.arange(1.5, 3.6, 0.1),
                ax=ax,
            )
    # delete both legends and remake my own
    for ax in axs:
        ax.get_legend().remove()
    # make my own legend
    axs[0].plot([], [], color=sns.color_palette("binary")[5], label="3 μm", linewidth=2)
    axs[0].plot([], [], color=sns.color_palette("binary")[2], label="13 μm", linewidth=2)
    # put legenbd to right of figure
    axs[0].axhspan(2.805, 3.005, color="red", alpha=0.2, label="cathode")
    axs[1].axhspan(2.805, 3.005, color="red", alpha=0.2, label="_")
    axs[0].axhspan(2.005, 2.205, color="blue", alpha=0.2, label="anode")
    axs[1].axhspan(2.005, 2.205, color="blue", alpha=0.2, label="_")

    axs[0].legend(loc="center left", bbox_to_anchor=(2.2, 0.5))
    axs[0].set_xlim(reversed(axs[0].get_xlim()))
    axs[0].set_ylabel("Activation Location (cm)\n(at threshold)")
    axs[0].set_title("extrusion", pad=-10)
    axs[1].set_title("true-3D", pad=-10)
    plt.suptitle(f'{waveform}', y=1)
    plt.gcf().set_size_inches(2, 2)
    plt.subplots_adjust(wspace=0.1)
# %% activation position
for waveform in newdefdat_allwf.waveform.unique():
    print(waveform)
    plt.figure()
    sns.set(font_scale=1.75, style="white")
    newthreshz = newdefdat_allwf.query("deformation=='Structural'")
    newthreshz["activation_zpos"] = newthreshz["activation_zpos"] / 10000
    fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
    for nerve_label in pd.unique(newthreshz.nerve_label):
        for ax, modeltype in zip(axs, ["extrusion", "true-3D"]):
            for fiber_diam in newthreshz.fiber_diam.unique():
                color = {3: 'black', 13: 'red'}[fiber_diam]
                lw = {3: 3, 13: 1.5}[fiber_diam]
                g = sns.histplot(
                    data=newthreshz.query(
                        f"fiber_diam ==@fiber_diam and waveform==@waveform and contact in {main_comparison} and nerve_label=='{nerve_label}' and type=='{modeltype}'"
                    ).rename(columns={"nerve_label": "Sample"}),
                    y="activation_zpos",
                    # hue="fiber_diam",
                    # hue='Sample',
                    # facet_kws={'sharex': False},
                    # kind='kde',
                    common_norm=False,
                    # legend=False,
                    # multiple="fill",
                    color=color,
                    lw=lw,
                    # color='black',alpha=0.5,
                    element="poly",
                    fill=False,
                    bins=np.arange(1.5, 3.6, 0.1),
                    ax=ax,
                )
    # # delete both legends and remake my own
    # for ax in axs:
    #     ax.get_legend().remove()
    # make my own legend
    axs[0].plot([], [], color='k', label="3 μm", linewidth=1)
    axs[0].plot([], [], color='k', label="13 μm", linewidth=3)
    # put legenbd to right of figure
    axs[0].axhspan(2.805, 3.005, color="red", alpha=0.2, label="cathode")
    axs[1].axhspan(2.805, 3.005, color="red", alpha=0.2, label="_")
    axs[0].axhspan(2.005, 2.205, color="blue", alpha=0.2, label="anode")
    axs[1].axhspan(2.005, 2.205, color="blue", alpha=0.2, label="_")

    axs[0].legend(loc="center left", bbox_to_anchor=(2.2, 0.5))
    axs[0].set_xlim(reversed(axs[0].get_xlim()))
    axs[0].set_ylabel("Activation Location (cm)\n(at threshold)")
    axs[0].set_title("extrusion")
    axs[1].set_title("true-3D")

# %% Thresholds plot
newdefdat_allwf.query(f'deformation == "Structural" and contact in {cath_comparison}')
# Group by fiber diameter, nerve label, and type
group_cols = ['fiber_diam', 'nerve_label', 'type']
# Calculate the median threshold for Asymmetric & Bipolar waveform
median_asym_bipolar = (
    newdefdat_allwf.query(
        'deformation == "Structural" and contact in @cath_comparison and waveform == "Asymmetric & Bipolar"'
    )
    .groupby(group_cols)['threshold']
    .mean()
    .reset_index()
    .rename(columns={'threshold': 'median_threshold'})
)
# Merge with the original data
normalized_data = (
    newdefdat_allwf.query('deformation == "Structural" and contact in @cath_comparison')
    .reset_index()
    .merge(median_asym_bipolar, on=group_cols, how='left')
)
# Normalize the threshold by the median of Asymmetric & Bipolar
normalized_data['threshold_normalized'] = normalized_data['threshold'] / normalized_data['median_threshold']
g = sns.catplot(
    kind='point',
    data=normalized_data,
    y='threshold_normalized',
    x='waveform',
    hue='fiber_diam',
    row='type',
    estimator='mean',
    dodge=True,
    # errorbar=('pi',100)
    # palette=wfpal,
    # sharey=False,
)
# g.set_titles(col_template='D: {col_name} μm')
# g.set_xlabels('Nerve')
# g.set_ylabels('Normalized Threshold')
# sns.move_legend(g,[0.75,0.4])
for ax in g.axes.ravel():
    plt.sca(ax)
    plt.axhline(1, color='black', ls='--')
    plt.xticks(rotation=20)
plt.gcf().set_size_inches([4, 4])
plt.subplots_adjust(hspace=0.2)
import matplotlib.pyplot as plt
import numpy as np

# %%
import pandas as pd
import seaborn as sns

# Define a "fiber key" that identifies each fiber (excluding waveform).
fiber_key = [
    "contact",
    "fiber_diam",
    "nerve_label",
    "type",
    "deformation",
    "master_fiber_index",
    # add others if needed, but NOT "waveform"
]

# 1) Extract exactly one baseline threshold (waveform="Asymmetric & Bipolar")
df_baseline = (
    newdefdat_allwf.query('waveform == "Asymmetric & Bipolar"')
    .drop_duplicates(subset=fiber_key)  # ensure exactly 1 row per fiber_key
    .copy()
)

# Rename the threshold column to "baseline_threshold"
df_baseline.rename(columns={"threshold": "baseline_threshold"}, inplace=True)

# 2) Merge that baseline with the full data
df_merged = pd.merge(newdefdat_allwf, df_baseline[fiber_key + ["baseline_threshold"]], on=fiber_key, how="left")

# 3) Compute normalized threshold by dividing by the single baseline
df_merged["threshold_normalized"] = df_merged["threshold"] / df_merged["baseline_threshold"]

# Example: Filter down if needed (e.g., only "Structural" + certain contacts)
df_plot = df_merged.query('deformation == "Structural" and contact in @cath_comparison')

# 4) Plot the normalized threshold vs. waveform, etc.
g = sns.catplot(
    kind='point',
    data=df_plot,
    x='waveform',
    y='threshold_normalized',
    hue='fiber_diam',  # or whatever makes sense
    row='type',
    estimator='mean',  # or median
    dodge=True,
)
# Optionally draw a horizontal line at 1, since that's the baseline
for ax_row in g.axes:
    for ax in ax_row:
        ax.axhline(1, ls='--', color='gray')
        plt.sca(ax)
        plt.xticks(rotation=20)

plt.gcf().set_size_inches(5, 6)
g.set_ylabels('Threshold change')
plt.tight_layout()
plt.show()


# %% Thresholds plot
newdefdat_allwf.query(f'deformation == "Structural" and contact in {cath_comparison}')
g = sns.catplot(
    data=newdefdat_allwf.query(f'deformation == "Structural" and contact in {cath_comparison}').reset_index(),
    y='threshold',
    hue='waveform',
    kind='box',
    x='nerve_label',
    col='fiber_diam',
    row='type',
    # col_wrap=2,
    height=4,
    aspect=1,
    palette=wfpal,
    sharey=False,
    # legend=False
)
g.set_titles(col_template='D: {col_name} μm')
g.set_xlabels('Nerve')
g.set_ylabels('Threshold (mA)')
sns.move_legend(g, [0.75, 0.4])
plt.gcf().set_size_inches([8, 6])
