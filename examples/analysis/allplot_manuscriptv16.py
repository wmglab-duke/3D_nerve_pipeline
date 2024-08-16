"""Created on Wed Mar  6 12:35:05 2024.

@author: dpm42
"""

import json
import os

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D
from scipy import stats
from scipy.stats import pearsonr, sem, variation

os.chdir("../../")
import sys

from src.core import Sample
from src.core.query import Query
from src.utils import Object

sys.path.pop(-2)
from src.core.plotter import datamatch_merge

mpl.rcParams["figure.dpi"] = 400
sns.set_style("whitegrid")
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
    with open(f"examples/analysis/{infile}.json") as f:
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


gogo = "initial"

cath_comparison = ["cathodic", "3D", 3]
center_comparison = ["center", "3D", 3]
an_comparison = ["anodic", "3D", 3]
comparisons = [cath_comparison, center_comparison, an_comparison]
main_comparison = cath_comparison  # NOTE: changed from center

simNUM = '333'

pal2d3d = ["#d95f02", "#7570b3"]

# code whole file to optionally run on 3 or 10 so that I can run all monpolar data
# %% set which comparison to run
print("threshload")
# base data
threshload = pd.read_csv(f"thresh_unmatched_sim{simNUM}_og.csv")

center = threshload["sample"].astype(str).str[2] == "1"
threshload.loc[center, "contact"] = "center"
threedcontact = threshload["sample"].astype(str).str[2] == "3"
threshload.loc[threedcontact, "contact"] = "3D"
# set all inners and outers where 'type' is 3D to 0
threshload.loc[(threshload["type"] == "3D"), "inner"] = 0
threshload.loc[(threshload["type"] == "3D"), "outer"] = 0
# %% inners need to match the cathodic leading contact
# Query the rows with type '3D'
df_3d = threshload.query("type == '3D'")

# Query the rows with type '2D' and sample ending in '1'
df_2d = threshload.query(f"type == '2D' and contact in {main_comparison}")

# Merge the 3D and 2D data, keeping track of original row indices
merged_df = pd.merge(
    df_3d,
    df_2d,
    on=["nerve_label", "master_fiber_index", "nsim"],
    suffixes=("_3d", "_2d"),
    how="left",
)  # TODO remove this how

# Update the 'inner', 'outer', and 'fiber' columns in the original DataFrame
threshload.loc[df_3d.index, "inner"] = merged_df["inner_2d"].values
threshload.loc[df_3d.index, "outer"] = merged_df["outer_2d"].values
threshload.loc[df_3d.index, "fiber"] = merged_df["fiber_2d"].values
# %%
if gogo == "initial":  # remove all where nerve_label length is >2
    threshload = threshload[threshload["nerve_label"].str.len() < 3]
threshload["type"] = threshload["type"].replace({"2D": "extrusion", "3D": "true-3D"})
threshload = addpwfd(threshload, str(simNUM), infile="plotconfig_og")
threshload["fiber_diam"] = threshload["fiber_diam"].astype(int)
# elif gogo=="deformedasc": #remove all where nerve_label does not contain "asc"
#     threshload = threshload[threshload['nerve_label'].str.contains("asc")]
#     #check the third digit of the sample number is 2, in that case, "contact" is cathodic. If 0, anodic
threshload["contact"] = (
    threshload["sample"].astype(str).str[2].replace({"0": "anodic", "2": "cathodic", "3": "3D", "1": "center"})
)
threshload.sort_values(by="sample")

# %%Setup
print("matched")
# true-3D and extrusion trhresholds as different col same row
matched = datamatch_merge(
    threshload.query('type=="extrusion"'),
    threshload.query('type=="true-3D"'),
    "threshold",
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns="type")
# for matched boxplots
repeated = matched.melt(
    id_vars=[x for x in matched.columns if "threshold" not in x],
    value_vars=["threshold3d", "threshold"],
    var_name="type",
    value_name="threshold",
)
repeated["type"] = repeated["type"].replace({"threshold": "extrusion", "threshold3d": "true-3D"})
repeated["type"] = pd.Categorical(repeated["type"], categories=["extrusion", "true-3D"], ordered=True)
# %%
print("dose-response")
# dose response
drdat = threshload.copy().rename(columns={"sample": "samplenum", "type": "modeltype"})
drdat = calculate_dose_response(
    drdat,
    "threshold",
    "percent_activated",
    grouping_columns=["modeltype", "samplenum", "fiber_diam", "sim"],
)
drdat.sort_values("modeltype", inplace=True)
drmatch = datamatch_merge(
    drdat.query('modeltype=="extrusion"'),
    drdat.query('modeltype=="true-3D"'),
    "percent_activated",
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns="modeltype")

# %% BEGIN DEFORMATION ANALYSIS
print("deformation")
# match fascicle thresholds
sns.set(font_scale=1.75)
sns.set_style("whitegrid")
threshes = pd.concat(
    [
        pd.read_csv(f"thresh_unmatched_sim{simNUM}_og.csv"),
        pd.read_csv(f"thresh_unmatched_sim{simNUM}_def.csv"),
    ]
)
# remove all rows where nerve label contains "asc" and "sample" contains "3"
threshes = threshes[~((threshes["nerve_label"].str.contains("asc")) & (threshes["sample"].astype(str).str[2] == "3"))]
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
newdefdat["nerve_label"] = pd.Categorical(newdefdat["nerve_label"], categories=["2L", "3R", "5R", "6R"], ordered=True)
newdefdat["deformed"] = newdefdat["deformation"] != "Undeformed"
deftomatch = newdefdat.copy()

# remove unused colors from palette
defpal = [sns.color_palette("colorblind")[ind] for ind in [0, 2, 3, 5]]
defdefcomp = newdefdat.query('type=="true-3D" and deformation != "ASCENT"')
defdefcomp["deformed"] = defdefcomp["deformation"] != "Undeformed"
defsamples = ["2L", "3R", "5R", "6R"]
# %% inners need to match the cathodic leading contact
definnermatch = deftomatch.query('deformation!="ASCENT"').reset_index(drop=True)

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
    thismatch = deftomatch.query("deformation==@deftype")
    matchednow = datamatch_merge(
        thismatch.query('type=="extrusion"'),
        thismatch.query('type=="true-3D"'),
        "threshold",
        merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
    ).drop(columns="type")
    matchednow["deformation"] = deftype
    concats.append(matchednow)
concats = pd.concat(concats)
concats["deformation"] = pd.Categorical(concats["deformation"], categories=["Undeformed", "Structural"], ordered=True)

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
    deftomatch.query('type=="extrusion" and deformation!="ASCENT"'),
    deftomatch.query('type=="true-3D" and deformation!="ASCENT"'),
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
sys.exit("prepdone")
# %% dose-response
sns.set(context='paper', style='ticks')
alldr = defdr.query(f"nerve_label in {defsamples}").query("nerve_label=='6R'")
alldr["deformed"] = alldr["deformation"] != "Undeformed"
allstore = []
for contact in alldr.contact.unique():
    if contact == "3D":
        continue
    thisdat = alldr.query('contact in [@contact, "3D"]')
    thisdat["contact"] = contact
    allstore.append(thisdat)
alldr = pd.concat(allstore)
alldr["deformation"] = pd.Categorical(
    alldr["deformation"],
    ordered=True,
    categories=["Undeformed", "Structural", "ASCENT"],
)

plt.figure()
g = sns.relplot(
    kind="line",
    style="deformation",
    data=alldr.query("fiber_diam in [3]"),
    y="percent_activated",
    x="threshold",
    units="nerve_label",
    hue="modeltype",
    palette=pal2d3d,
    estimator=None,
    linewidth=2,
    facet_kws={"sharex": True, "margin_titles": True},
    row="deformed",
    col="contact",
)

g.legend.set_title("")

for ax in g.axes.ravel():
    ax.set_ylim([0, 1])
    ax.set_xlim([0, None])
g.set_ylabels("Proportion fibers active")
g.set_titles(col_template="{col_name}", row_template='')
plt.subplots_adjust(hspace=0.15)

g.set_xlabels("Threshold (mA)")
# change the line width for the legend
for line, l in zip(g.legend.get_lines(), g.legend.get_texts()):
    line.set_linewidth(2.0)
    # if l.get_text() in ['deformation', 'modeltype']:
    #     l.set_text('')
for ax in g.axes.ravel():
    for loc in [0.1, 0.5, 0.9]:
        ax.axhline(loc, color="gray", linestyle="-", alpha=0.5, linewidth=1)
sns.move_legend(g, [0.75, 0.2], facecolor='white', framealpha=1, frameon=True, edgecolor='black')
g.axes[0][0].set_ylabel("Active fibers (%)\nUndeformed")
g.axes[1][0].set_ylabel("Active fibers (%)\nDeformed")
# rownames(g, row_template="{row_name}-slice\nProportion fibers active")
g.axes[0][0].set_yticks([0, 0.1, 0.5, 0.9, 1], ['0', '10', '50', '90', '100'])
g.fig.set_size_inches(4, 2.5)
# %% threshold vals
for stringdat in ["Undeformed", "Structural"]:
    # deformstate = deformation!="Undeformed"
    thisdat = datamatch_merge(
        newdefdat.query(f'type=="extrusion" and deformation=="{stringdat}" and contact in {main_comparison}'),
        newdefdat.query(f'type=="true-3D" and deformation=="{stringdat}" and contact in {main_comparison}'),
        "threshold",
        merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
    ).drop(columns="type")
    for diam in [3, 13]:
        dadata = thisdat.query(f"fiber_diam==@diam and contact in {main_comparison}")
        perc = sum(dadata.threshold > dadata.threshold3d) / len(dadata.threshold)
        print(f"Percent higher thresh for {diam}:", round(perc, 3))
        res = np.abs(dadata.threshold3d - dadata.threshold)
        print(
            f"Mean abs residual for {diam}:",
            round(np.mean(res), 3),
            "sem:",
            round(sem(res), 3),
        )
        print(
            f"Max,min,mean residual for {diam}:",
            round(np.max(res), 3),
            round(np.min(res), 3),
            round(np.mean(res), 3),
            "+-",
            round(sem(res), 3),
        )
        mean_simthresh = np.mean(np.concatenate([dadata.threshold, dadata.threshold3d]))
        print(f"Mean sim thresh for {diam}:", round(mean_simthresh, 3))
        # max, min, mean, sem as a percentage of mean_simthresh
        print(
            f"Max,min,mean,sem as a percentage of mean_simthresh for {diam}:",
            round(np.max(res) / mean_simthresh, 3),
            round(np.min(res) / mean_simthresh, 3),
            round(np.mean(res) / mean_simthresh, 3),
            "+-",
            round(sem(res) / mean_simthresh, 3),
        )
        res = dadata.threshold3d - dadata.threshold
print(
    f"Difference between means for {diam}:",
    round(np.mean(res), 3),
    "sem:",
    round(sem(res), 3),
)
# %% dose-response onset sat deformed
sns.set(context='paper', style='whitegrid')
thisheredata = []
for comparison in [an_comparison, center_comparison, cath_comparison]:
    subdat = newdefdat.query(f"contact in {comparison}")
    plt.figure()
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
            "deformation",
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
    compiled_data["contact"] = comparison[0]
    thisheredata.append(compiled_data)

allcomp = pd.concat(thisheredata)
compmatch = datamatch_merge(
    allcomp.query('type=="extrusion"'),
    allcomp.query('type=="true-3D"'),
    "threshold",
    merge_cols=[
        "model",
        "sim",
        "nerve_label",
        "nsim",
        "level",
        "deformation",
        "contact",
    ],
).drop(columns="type")

g = sns.relplot(
    data=compmatch.query("fiber_diam in [3,13]").rename(columns={"nerve_label": "Sample"}),
    kind="scatter",
    row="contact",
    col="fiber_diam",
    style="deformation",
    x="threshold",
    y="threshold3d",
    hue="level",
    s=60,
    palette="gnuplot2",
    facet_kws={"sharex": "col", "sharey": "col", "margin_titles": True},
    edgecolor="black",
    linewidth=1,
    alpha=1,
)
# TODO Clean up this calc
for ax in g.axes.ravel():
    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim], [0, lim], "--k", linewidth=2, label="unity line", zorder=-1)  # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])

    # ax.set_yticks(ax.get_xticks())
g.set_titles(col_template="D: {col_name} μm", row_template="")
g.set_xlabels("extrusion Threshold (mA)")
g.set_ylabels("true-3D Threshold (mA)")
# plt.legend()
compmatch["threshdiff"] = compmatch.threshold3d - compmatch.threshold
compmatch.deformation = pd.Categorical(compmatch.deformation, categories=["Structural", "ASCENT"], ordered=True)
g = sns.relplot(
    data=compmatch.query('deformation!= "Undeformed" and contact=="cathodic" and fiber_diam in [3,13]').rename(
        columns={"nerve_label": "Sample"}
    ),
    kind="scatter",
    # row="contact",
    row="fiber_diam",
    style="deformation",
    x="threshold",
    y="threshold3d",
    hue="level",
    # s=60,
    # size='level',
    palette="gnuplot2",
    facet_kws={"sharex": False, "sharey": False, "margin_titles": True},
    edgecolor="black",
    # linewidth=1,
    alpha=1,
    # col_wrap=2,
)
# TODO Clean up this calc
for ax in g.axes.ravel():
    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim], [0, lim], "--k", linewidth=2, label="unity line", zorder=-1)  # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])

    # ax.set_yticks(ax.get_xticks())
plt.subplots_adjust(hspace=0.2)
sns.move_legend(g, [0.68, 0.3])
g.set_titles(col_template="", row_template="")
g.set_xlabels("extrusion threshold (mA)")
# g.set_ylabels("true-3D Threshold (mA)")
rownames(g, "true-3D threshold (mA)\nD: {row_name} μm")
plt.gcf().set_size_inches(2.5, 3.5)
# %% dose-response onset sat deformed
peses = []
pemeans = []
onsets_sats = {}
for comparison in comparisons:
    for stringdat in ["Undeformed", "Structural", "ASCENT"]:
        thiscontact = comparison[0]
        subdat = newdefdat.query(f"deformation=='{stringdat}' and contact in {comparison}")
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
                    pe_res = pe_noabs(b[0], a[0], doabs=True)
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
        assert thiscontact != np.nan
        peses.append(pes)
        print("Max 3 um", stringdat, np.amax(pes.query("fiber_diam==3").pe))
        print("Max 13 um", stringdat, np.amax(pes.query("fiber_diam==13").pe))

        # now calculate percent error for population onset and saturation
        pemean = []
        for level in ["onset", "half", "saturation"]:
            for fiber_diam in compiled_data.fiber_diam.unique():
                a = compiled_data.query(f"level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam}")[
                    "threshold"
                ].values
                b = compiled_data.query(f"level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}")[
                    "threshold"
                ].values
                pe_res = pe_noabs(np.median(b), np.median(a), doabs=True)
                pemean.append({"level": level, "fiber_diam": fiber_diam, "pe": pe_res})

        pemean = pd.DataFrame(pemean)
        pemean["deformation"] = stringdat
        pemean["contact"] = thiscontact
        pemeans.append(pemean)
        print("Max 3 um median", stringdat, np.amax(pemean.query("fiber_diam==3").pe))
        print("Max 13 um median", stringdat, np.amax(pemean.query("fiber_diam==13").pe))

        onsets_sats[stringdat] = compiled_data

sns.set(style="whitegrid", context="paper")
allpes = pd.concat(peses)
# allpes['level'].replace({
#     'onset':'10','half':'50','saturation':'90'},inplace=True)
# allpes['deformation'].replace({
#     'ASCENT':'ASCENT','Structural':'Structural','Undeformed':'Undef.'},inplace=True)
# allpes['comb'] = allpes['level'].astype(str) + '\n' + allpes['deformation']
# allpes.sort_values(by=['deformation','level'])
allpes["level"] = pd.Categorical(allpes.level, categories=["onset", "half", "saturation"], ordered=True)
allpes["deformation"] = pd.Categorical(
    allpes.deformation, categories=["Undeformed", "Structural", "ASCENT"], ordered=True
)
allpes["contact"] = pd.Categorical(allpes.contact, categories=["anodic", "center", "cathodic"], ordered=True)
allpes.sort_values(by=["contact", "deformation"], inplace=True)
allpes["comb"] = allpes["deformation"].astype(str) + "\n" + allpes["contact"].astype(str)

sns.set(context='paper', style='white', font_scale=1)
# allpes=allpes.query('deformation!="ASCENT"')
# allpes["deformation"] = allpes.deformation.replace({'Structural':'Deformed'})
# allpes["deformation"] = pd.Categorical(
#     allpes.deformation, categories=["Undeformed", "Deformed"], ordered=True
# )
# rewrite the above with map dataframe
allpes_less = allpes.copy()
allpes_less["deformation"] = pd.Categorical(allpes_less.deformation, categories=["Structural", "ASCENT"], ordered=True)
allpes_less["contact"] = pd.Categorical(allpes_less.contact, categories=["cathodic"], ordered=True)

g = sns.FacetGrid(
    data=allpes_less.query("fiber_diam==3"),
    palette=defpal,
    col="contact",
    row="deformation",
    # row_order=["saturation", "half", "onset"],
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x="level",
    y="pe",
    marker="o",
    alpha=0.4,
    units="nerve_label",
    estimator=None,
    palette=defpal,
    color='gray',
    # style='deformation'
)
g.map_dataframe(
    sns.lineplot, x="level", y="pe", marker="s", estimator="median", errorbar=None, color="black", label='median'
)
plt.gcf().set_size_inches(4, 4)
g.axes.ravel()[0].legend()
g.set_titles(col_template="{col_name}", row_template="")
rownames(g, row_template="{row_name}")
g.axes[1, 0].set_ylabel(f'Absolute Percent Difference (%)\n{g.axes[1,0].get_ylabel()}')
for ax in g.axes.ravel():
    plt.sca(ax)
    # plt.xticks(rotation=20)
    plt.xlabel("% active")
    ax.set_xticklabels(['10\n(onset)', '50\n(half)', '90\n(sat.)'])
plt.subplots_adjust(hspace=0.1, wspace=0)
# plt.gcf().set_size_inches(9,9)
plt.ylim(0, 80)
plt.xlim(-0.5, 2.5)
# rewrite the above with map dataframe
g = sns.FacetGrid(
    data=allpes_less,
    palette=defpal,
    col="deformation",
    row="contact",
    # row_order=["saturation", "half", "onset"],
    margin_titles=True,
)
g.map_dataframe(
    sns.barplot,
    x="level",
    y="pe",
    # marker="s",
    estimator="median",
    errorbar=None,
    hue="fiber_diam",
    palette=rdup,
    # markersize=3,
    # linewidth=1,
    # dodge=0.6,
)
g.map_dataframe(
    sns.stripplot,
    x="level",
    y="pe",
    # marker="s",
    # estimator="median",
    # errorbar=None,
    hue="fiber_diam",
    # palette=rdup,
    # markersize=3,
    linewidth=0.5,
    dodge=True,
    edgecolor='w',
    color='k',
    jitter=False,
    s=3,
)
# plt.gcf().set_size_inches(6, 6)
plt.legend(ncol=1, bbox_to_anchor=[1, 1], title='D (μm)', framealpha=0)
g.set_titles(col_template="{col_name}", row_template="")
# rownames(g, row_template="{row_name}")
g.set_ylabels('Absolute\nPercent Difference (%)')
# g.axes[1,0].set_ylabel(f'Absolute Percent Difference (%)\n{g.axes[1,0].get_ylabel()}')
for ax in g.axes.ravel():
    plt.sca(ax)
    # plt.xticks(rotation=20)
    plt.xlabel("% active")
    ax.set_xticklabels(['10\n(onset)', '50\n(half)', '90\n(sat.)'])
for ax in g.axes.ravel():
    ax.tick_params(pad=0)
plt.subplots_adjust(hspace=0.1, wspace=0)
# plt.ylim(0, 20)
plt.xlim(-0.5, 2.5)
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title="", handles=handles[:6], labels=labs[:6], bbox_to_anchor=[1, 1])
plt.gcf().set_size_inches(10 / 3, 5 / 3)
print(
    allpes.query(f'contact in {cath_comparison} and deformation=="Structural"')
    .groupby(['fiber_diam', 'level'])
    .median()
    .groupby(['level'])
    .max()
)
# %% dose-response onset sat deformed
peses = []
pemeans = []
onsets_sats = {}
for comparison in comparisons:
    for stringdat in ["Undeformed", "Structural", "ASCENT"]:
        thiscontact = comparison[0]
        subdat = newdefdat.query(f"deformation=='{stringdat}' and contact in {comparison}")
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
                    pe_res = pe_noabs(b[0], a[0], doabs=True)
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
        assert thiscontact != np.nan
        peses.append(pes)
        print("Max 3 um", stringdat, np.amax(pes.query("fiber_diam==3").pe))
        print("Max 13 um", stringdat, np.amax(pes.query("fiber_diam==13").pe))

        # now calculate percent error for population onset and saturation
        pemean = []
        for level in ["onset", "half", "saturation"]:
            for fiber_diam in compiled_data.fiber_diam.unique():
                a = compiled_data.query(f"level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam}")[
                    "threshold"
                ].values
                b = compiled_data.query(f"level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}")[
                    "threshold"
                ].values
                pe_res = pe_noabs(np.median(b), np.median(a), doabs=True)
                pemean.append({"level": level, "fiber_diam": fiber_diam, "pe": pe_res})

        pemean = pd.DataFrame(pemean)
        pemean["deformation"] = stringdat
        pemean["contact"] = thiscontact
        pemeans.append(pemean)
        print("Max 3 um median", stringdat, np.amax(pemean.query("fiber_diam==3").pe))
        print("Max 13 um median", stringdat, np.amax(pemean.query("fiber_diam==13").pe))

        onsets_sats[stringdat] = compiled_data

sns.set(style="whitegrid", context="paper")
allpes = pd.concat(peses)
# allpes['level'].replace({
#     'onset':'10','half':'50','saturation':'90'},inplace=True)
# allpes['deformation'].replace({
#     'ASCENT':'ASCENT','Structural':'Structural','Undeformed':'Undef.'},inplace=True)
# allpes['comb'] = allpes['level'].astype(str) + '\n' + allpes['deformation']
# allpes.sort_values(by=['deformation','level'])
allpes["level"] = pd.Categorical(allpes.level, categories=["onset", "half", "saturation"], ordered=True)
allpes["deformation"] = pd.Categorical(
    allpes.deformation, categories=["Undeformed", "Structural", "ASCENT"], ordered=True
)
allpes["contact"] = pd.Categorical(allpes.contact, categories=["anodic", "center", "cathodic"], ordered=True)
allpes.sort_values(by=["contact", "deformation"], inplace=True)
allpes["comb"] = allpes["deformation"].astype(str) + "\n" + allpes["contact"].astype(str)

sns.set(context='paper', style='white', font_scale=1)
# allpes=allpes.query('deformation!="ASCENT"')
# allpes["deformation"] = allpes.deformation.replace({'Structural':'Deformed'})
# allpes["deformation"] = pd.Categorical(
#     allpes.deformation, categories=["Undeformed", "Deformed"], ordered=True
# )
# rewrite the above with map dataframe
g = sns.FacetGrid(
    data=allpes.query("fiber_diam==3"),
    palette=defpal,
    col="contact",
    row="deformation",
    # row_order=["saturation", "half", "onset"],
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x="level",
    y="pe",
    marker="o",
    alpha=0.4,
    units="nerve_label",
    estimator=None,
    palette=defpal,
    color='gray',
    # style='deformation'
)
g.map_dataframe(
    sns.lineplot, x="level", y="pe", marker="s", estimator="median", errorbar=None, color="black", label='median'
)
plt.gcf().set_size_inches(4, 4)
g.axes.ravel()[2].legend(frameon=False)
g.set_titles(col_template="{col_name}", row_template="")
rownames(g, row_template="{row_name}")
g.axes[1, 0].set_ylabel(f'Absolute Percent Difference (%)\n{g.axes[1,0].get_ylabel()}')
for ax in g.axes.ravel():
    ax.tick_params(pad=0)
    plt.sca(ax)
    # plt.xticks(rotation=20)
    plt.xlabel("% active")
    ax.set_xticklabels(['10\n(onset)', '50\n(half)', '90\n(sat.)'])
    ax.set_yticks([0, 25, 50, 75])
plt.subplots_adjust(hspace=0.1, wspace=0)
# plt.gcf().set_size_inches(9,9)
plt.ylim(0, 80)
plt.xlim(-0.5, 2.5)
# rewrite the above with map dataframe
g = sns.FacetGrid(
    data=allpes,
    palette=defpal,
    col="contact",
    row="deformation",
    # row_order=["saturation", "half", "onset"],
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x="level",
    y="pe",
    marker="s",
    estimator="median",
    errorbar=None,
    hue="fiber_diam",
    palette=rdup,
)
# plt.gcf().set_size_inches(6, 6)
plt.legend(ncol=3, bbox_to_anchor=[0.8, 1], title='D (μm)', framealpha=1)
g.set_titles(col_template="{col_name}", row_template="")
rownames(g, row_template="{row_name}")
g.axes[1, 0].set_ylabel(f'Absolute Percent Difference (%)\n{g.axes[1,0].get_ylabel()}')
for ax in g.axes.ravel():
    plt.sca(ax)
    # plt.xticks(rotation=20)
    plt.xlabel("% active")
    ax.set_xticklabels(['10\n(onset)', '50\n(half)', '90\n(sat.)'])
plt.subplots_adjust(hspace=0.1, wspace=0)
plt.ylim(0, 80)
plt.xlim(-0.5, 2.5)
plt.gcf().set_size_inches(4, 4)
# %%newfig dose-response
sns.set(style='white', context='paper')
g = sns.FacetGrid(
    data=allpes.query("fiber_diam in [3,13]"),
    palette=defpal,
    col="fiber_diam",
    row="deformation",
    # row_order=["saturation", "half", "onset"],
    margin_titles=True,
)
g.map_dataframe(
    sns.swarmplot,
    x="contact",
    y="pe",
    dodge=0.2,
    hue='level',
    s=2.5,
    legend=False,
    palette='gnuplot2',
    alpha=1,
    linewidth=0.5,
    edgecolor='w',
)
g.map_dataframe(
    sns.pointplot,
    x="contact",
    y="pe",
    marker="s",
    hue='level',
    estimator="median",
    errorbar=None,
    color="black",
    legend=True,
    palette='gnuplot2',
    dodge=0.5,
    markersize=3.5,
    linewidth=1,
)
plt.ylim(None, 80)
g.axes[0, 1].legend(frameon=False)
g.set_titles(col_template="D: {col_name} μm", row_template="")
rownames(g, row_template="{row_name}")
for ax in g.axes.ravel():
    ax.tick_params(pad=0)
    plt.sca(ax)
    # plt.xticks(rotation=20)
    # plt.xlabel("% active")
    plt.xticks(rotation=45)
g.axes[1, 0].set_ylabel(f'Absolute Percent Difference (%)\n{g.axes[1,0].get_ylabel()}')
plt.gcf().set_size_inches(3, 4)

# %%
g = sns.FacetGrid(
    data=allpes.query("fiber_diam==3"),
    palette=defpal,
    col="level",
    row="deformation",
    # row_order=["saturation", "half", "onset"],
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x="contact",
    y="pe",
    marker="o",
    alpha=0.4,
    units="nerve_label",
    estimator=None,
    palette=defpal,
    color='gray',
    # style='deformation'
)
g.map_dataframe(
    sns.lineplot, x="contact", y="pe", marker="s", estimator="median", errorbar=None, color="black", label='median'
)
plt.gcf().set_size_inches(4, 4)
g.axes.ravel()[0].legend()
g.set_titles(col_template="{col_name}", row_template="")
rownames(g, row_template="{row_name}")
g.axes[1, 0].set_ylabel(f'Absolute Percent Difference (%)\n{g.axes[1,0].get_ylabel()}')
for ax in g.axes.ravel():
    plt.sca(ax)
    # plt.xticks(rotation=20)
    plt.xlabel("% active")
    ax.set_xticklabels(['anodic', 'center', 'cathodic'], rotation=20)
for i, name in enumerate(['onset (10%)', 'half (50%)', 'saturation (90%)']):
    g.axes[0, i].set_title(name)
plt.subplots_adjust(hspace=0.1, wspace=0)
# plt.gcf().set_size_inches(9,9)
plt.ylim(None, 80)
plt.xlim(-0.5, 2.5)
# %%
g = sns.FacetGrid(
    data=allpes.query("fiber_diam==3"),
    palette=defpal,
    col="level",
    row="contact",
    # row_order=["saturation", "half", "onset"],
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x="deformation",
    y="pe",
    marker="o",
    alpha=0.4,
    units="nerve_label",
    estimator=None,
    palette=defpal,
    color='gray',
    # style='deformation'
)
g.map_dataframe(
    sns.lineplot, x="deformation", y="pe", marker="s", estimator="median", errorbar=None, color="black", label='median'
)
for ax in g.axes.ravel():
    plt.sca(ax)
    plt.xticks(rotation=20)
plt.gcf().set_size_inches(4, 4)
# smol
g = sns.FacetGrid(
    data=allpes.query("fiber_diam in [3,13] and contact in ['cathodic',3]"),
    palette=defpal,
    col='fiber_diam',
    # row_order=["saturation", "half", "onset"],
    margin_titles=True,
)
# g.map_dataframe(
#     sns.swarmplot, x="contact", y="pe",dodge=True,
#     hue='level',s=4,
# legend=False,
#     palette = 'gnuplot2',alpha=0.5
# )
g.map_dataframe(
    sns.barplot, x="deformation", y="pe", hue='level', estimator="median", errorbar=None, palette='gnuplot2'
)
plt.legend()
# plt.ylim(None, 80)
for ax in g.axes.ravel():
    plt.sca(ax)
    plt.xticks(rotation=20)
plt.gcf().set_size_inches(4, 3)
plt.gca().set_ylabel('Absolute Percent Difference (%)')
# %% now for only cathodic calculate the median and plot barplot with values
specdat = allpes.query(f'contact in {cath_comparison}')
specdat = specdat.groupby(['deformation', 'level', 'fiber_diam']).agg({'pe': 'median'}).reset_index()
g = sns.catplot(
    row='deformation', col='level', x='fiber_diam', y='pe', data=specdat, kind='bar', palette='RdPu', margin_titles=True
)
# with values
for ax in g.axes.ravel():
    for i in ax.containers:
        ax.bar_label(i, fmt='%.1f')
plt.gcf().set_size_inches(4, 4)
# %% change in threshold after deformation
diffs = []
# calculate the change in onset and saturation for each sample (compared between undeformed and Structural deformed)
for nerve_label in pd.unique(onsets_sats["Undeformed"].nerve_label):
    for fiber_diam in pd.unique(onsets_sats["Undeformed"].fiber_diam):
        for level in ["onset", "saturation"]:
            undeformed = (
                onsets_sats["Undeformed"]
                .query(
                    f"type=='true-3D' and nerve_label == '{nerve_label}' and fiber_diam == {fiber_diam} and level == '{level}'"
                )["threshold"]
                .values[0]
            )
            deformed = (
                onsets_sats["Structural"]
                .query(
                    f"type=='true-3D' and nerve_label == '{nerve_label}' and fiber_diam == {fiber_diam} and level == '{level}'"
                )["threshold"]
                .values[0]
            )
            # print(f"{nerve_label} {fiber_diam} {level} {deformed-undeformed/undeformed}")
            diffs.append(
                {
                    "nerve_label": nerve_label,
                    "fiber_diam": fiber_diam,
                    "level": level,
                    "diff": 100 * (deformed - undeformed) / undeformed,
                }
            )
diffs = pd.DataFrame(diffs)
# print mean, std, min, max
print("Mean difference", np.mean(diffs["diff"]))
print("Std difference", np.std(diffs["diff"]))
print("Min difference", np.min(diffs["diff"]))
print("Max difference", np.max(diffs["diff"]))
# %% activation order calc
newdefdr = defdr.copy()
newdefdr = newdefdr.query("deformation=='Structural'").sort_values("modeltype")
sns.set(font_scale=1.25)
sns.set_style("whitegrid")
newdefdr["percent_activated2d"] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newdefdr.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        'modeltype == "true-3D" and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim and master_fiber_index == @row.master_fiber_index and deformation==@row.deformation'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert val is not np.nan
    newdefdr.loc[row.Index, "percent_activated3d"] = val
newdefdr["contact"] = newdefdr["contact"].replace(
    {
        "3D": "true-3D",
        "cathodic": "cathodic-leading\nextrusion",
        "anodic": "anodic-leading\nextrusion",
        "center": "center\nextrusion",
    }
)
newdefdr["contact"] = pd.Categorical(
    newdefdr["contact"],
    categories=[
        "true-3D",
        "anodic-leading\nextrusion",
        "center\nextrusion",
        "cathodic-leading\nextrusion",
    ],
    ordered=True,
)
# %% plot activation order
sns.set(font_scale=1, context='paper', style="white")
plotactidata = newdefdr.query("fiber_diam in [3] and nerve_label=='3R'")
plt.figure()
g = sns.catplot(
    kind="swarm",
    row="fiber_diam",
    data=plotactidata,
    y="percent_activated",
    x="contact",
    units="nerve_label",
    palette="plasma",
    hue="percent_activated3d",
    # hue='inner',
    estimator=None,
    linewidth=0,
    # facet_kws={"margin_titles": True},
    s=10,
)
plt.subplots_adjust(top=0.87)
# plt.suptitle(stringdat, x=0.37)
g.set_titles(row_template="", col_template="{col_name}")
g.axes[0][0].set_xlabel("")
g.axes[0][0].set_ylabel("Proportion fibers activated")
g.set_xlabels("")
g.legend.remove()
norm = plt.Normalize(0, 1)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])
plt.axvline(0.5, linestyle="--", color="black", alpha=0.5)
# Remove the legend and add a colorbar
g.figure.colorbar(
    sm,
    ax=g.axes.ravel().tolist(),
    aspect=10,
    shrink=0.8,
    label="Proportion true-3D\nfibers activated",
    pad=0.1,
).ax.yaxis.set_ticks_position("left")
g.fig.set_size_inches(14, 4)
for i, con in enumerate(newdefdr.contact.sort_values().unique()[1:]):
    shortdat = plotactidata.query("contact==@con")
    data2d = shortdat.sort_values("percent_activated").master_fiber_index
    data3d = shortdat.sort_values("percent_activated3d").master_fiber_index
    rc = compute_reorder_cost(list(data2d), list(data3d))
    g.axes[0][0].text(i + 0.6, 1.065, f'AR: {round(rc, 3)}')
plt.gcf().set_size_inches(6, 2)
# %%plot FASR
thiscomp = cath_comparison


def recruitment_cost_more(data):
    """Calculate recruitment cost for each inner.

    Recruitment cost is defined as the ratio of number of stimulated
    off-target fibers to total number of off-target fibers. From
    https://iopscience.iop.org/article/10.1088/1741-2560/10/3/036010
    """
    for inner in pd.unique(data["inner"]):
        # get threshold for inner
        inner_data = data.query(f"inner == {inner}")["threshold"]
        assert len(data) > 200 & len(data) < 250
        inner_thresh = np.amax(inner_data)
        # get all off-target fiber thresholds
        off_thresh = data.query(f"inner != '{inner}'")["threshold"]
        # calculate recruitment cost
        cost = np.sum(off_thresh <= inner_thresh) / len(off_thresh)
        data.loc[data["inner"] == inner, "RC"] = cost
        # fascicle selectivity ratio is 1-RC
        data.loc[data["inner"] == inner, "FASR"] = 1 - cost
        fasr_dict = {
            "active_src_index": data["active_src_index"].iloc[0],
            "fiber_diam": data["fiber_diam"].iloc[0],
            "type": data["type"].iloc[0],
            "inner": inner,
            "RC": cost,
            "FASR": 1 - cost,
            "nerve_label": data["nerve_label"].iloc[0],
            "deformation": data["deformation"].iloc[0],
            "contact": data["contact"].iloc[0],
        }
        yield fasr_dict


sns.set(font_scale=1.5)
sns.set_style("whitegrid")
FASRS = []
for deformation in definnermatch.deformation.unique():
    for nerve in pd.unique(definnermatch["nerve_label"]):
        for n in definnermatch.fiber_diam.unique():
            fasrdatto = definnermatch.query(f"contact in {thiscomp} and deformation==@deformation").query(
                f'nerve_label=="{nerve}" and fiber_diam=={n}'
            )
            for typ in fasrdatto.type.unique():
                FASRS.extend(recruitment_cost_more(fasrdatto.query("type==@typ")))

allfasr = pd.DataFrame(FASRS)

g = sns.catplot(
    data=allfasr.query('deformation=="Structural"'),
    x="type",
    y="FASR",
    hue="fiber_diam",
    dodge=True,
    palette="RdPu",
    kind="box",
)
g.set_xlabels("")
g.legend.set_title("D: μm")

# %%
fasnew = datamatch_merge(
    allfasr.query('type=="extrusion"'),
    allfasr.query('type=="true-3D"'),
    "FASR",
    merge_cols=["inner", "nerve_label", "deformation", "fiber_diam"],
).drop(columns="type")
fasnew["FASR-diff"] = fasnew["FASR3d"] - fasnew["FASR"]
sns.violinplot(data=fasnew, y="FASR-diff", x="fiber_diam", fill=True, linecolor="k", color="w")
sns.swarmplot(data=fasnew, y="FASR-diff", x="fiber_diam", color="k", s=2)
plt.ylim(-1, 1)
plt.xlabel("D (μm)")
plt.ylabel("FASR difference\n(true-3D minus extrusion)")

# %% calc stats
datahere = deftomatch.query('deformation in ["Structural","Undeformed"]')
mathere = concats.copy()
scores = []
FASRS = []
for comp in [center_comparison, cath_comparison, an_comparison]:
    for deformation in datahere.deformation.unique():
        threshdat = datahere.query(f"contact in {comp} and deformation==@deformation")
        for nerve in pd.unique(threshdat["nerve_label"]):
            for n in [3, 13]:
                shortdat = threshdat.query(f'nerve_label=="{nerve}" and fiber_diam=={n}')
                data2d = shortdat.query('type=="extrusion"').sort_values("threshold").master_fiber_index
                data3d = shortdat.query('type=="true-3D"').sort_values("threshold").master_fiber_index
                rc = compute_reorder_cost(list(data2d), list(data3d))
                concathere = mathere.query(
                    f'nerve_label=="{nerve}" and fiber_diam=={n}and contact in {comp} and deformation==@deformation'
                )
                data2d = concathere.threshold
                data3d = concathere.threshold3d
                ccc = concordance_correlation_coefficient(list(data2d), list(data3d))
                scores.append(
                    {
                        "sample": nerve,
                        "fiber_diam": n,
                        "score2d3d": rc,
                        "deformation": deformation,
                        "slice": comp[0],
                        "CCC": ccc,
                    }
                )
scoredat = pd.DataFrame(scores)
# %%
sns.set(context='paper', style='whitegrid')
g = sns.FacetGrid(data=scoredat, margin_titles=True, col="fiber_diam")
g.map_dataframe(
    sns.barplot,
    y="CCC",
    x="slice",
    hue="deformation",
    errorbar=None,
    palette="binary",
    order=['anodic', 'center', 'cathodic'],
    estimator='median',
)
g.map_dataframe(
    sns.stripplot,
    y="CCC",
    x="slice",
    hue="deformation",
    dodge=True,
    order=['anodic', 'center', 'cathodic'],
    color='black',
    edgecolor='white',
    linewidth=0.5,
    s=4,
)

# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel("CCC")
g.set_titles(col_template="D: {col_name} μm")
g.set_ylabels("CCC")
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
# plt.ylim(0, 0.6)
print(scoredat.groupby(["deformation"]).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
# g.fig.set_size_inches(8,4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel("extrusion slice")
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title="", handles=handles[:2], labels=labs[:2], bbox_to_anchor=[0.8, 1.5], ncol=2)
plt.gcf().set_size_inches(3, 1.5)
plt.ylim(0, 1)
# grab median vals by fiber diameter and print for cathodic structural deformation
print(scoredat.query('deformation=="Structural" and slice=="cathodic"').groupby(["fiber_diam"]).median())
# %%
sns.set(context='paper', style='whitegrid')
g = sns.FacetGrid(data=scoredat, margin_titles=True, col="fiber_diam")
g.map_dataframe(
    sns.barplot,
    y="score2d3d",
    x="slice",
    hue="deformation",
    errorbar=None,
    palette="binary",
    order=['anodic', 'center', 'cathodic'],
    estimator='median',
)
g.map_dataframe(
    sns.stripplot,
    y="score2d3d",
    x="slice",
    hue="deformation",
    dodge=True,
    order=['anodic', 'center', 'cathodic'],
    color='black',
    edgecolor='white',
    linewidth=0.5,
    s=4,
)

# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel("Activation Reordering")
g.set_titles(col_template="D: {col_name} μm")
g.set_ylabels("Activation Reordering")
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
# plt.ylim(0, 0.6)
print(scoredat.groupby(["deformation"]).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
# g.fig.set_size_inches(8,4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel("extrusion slice")
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title="", handles=handles[:2], labels=labs[:2], bbox_to_anchor=[0.8, 1.5], ncol=2)
plt.gcf().set_size_inches(3, 1.5)
plt.ylim(0, 1)
# add a horizontal line at y=0.45 for "completely random" and label with text in red
for ax in g.axes.ravel():
    ax.axhline(0.45, color='red', linestyle='--')
    ax.text(1.4, 0.47, 'random', color='red')
# print max and min values for cathodic structural deformation
print(
    scoredat.query('deformation=="Structural" and slice=="cathodic"')
    .groupby(["fiber_diam"])[["score2d3d"]]
    .agg([np.min, np.max])
)
# %% vals
g = sns.barplot(
    data=scoredat.query('slice=="cathodic"'), y='score2d3d', x="fiber_diam", estimator='median', hue='deformation'
)
for i in g.containers:
    g.bar_label(
        i,
    )

# %% remake the above with individual calls of histplot
sns.set(font_scale=1.75, style="whitegrid")
newthreshz = newdefdat.query("deformation=='Structural'")
newthreshz["activation_zpos"] = newthreshz["activation_zpos"] / 10000
fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
for nerve_label in pd.unique(newthreshz.nerve_label):
    for ax, modeltype in zip(axs, ["extrusion", "true-3D"]):
        g = sns.histplot(
            data=newthreshz.query(
                f"fiber_diam in [3,13] and contact in {main_comparison} and nerve_label=='{nerve_label}' and type=='{modeltype}'"
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
axs[0].axhspan(2.005, 2.205, color="blue", alpha=0.2, label="anode")
axs[1].axhspan(2.005, 2.205, color="blue", alpha=0.2, label="_")
axs[0].axhspan(2.805, 3.005, color="red", alpha=0.2, label="cathode")
axs[1].axhspan(2.805, 3.005, color="red", alpha=0.2, label="_")

axs[0].legend(loc="center left", bbox_to_anchor=(2.2, 0.5))
axs[0].set_xlim(reversed(axs[0].get_xlim()))
axs[0].set_ylabel("Activation Location (cm)")
axs[0].set_title("extrusion-100%")
axs[1].set_title("true-3D-100%")

# %% Percent Error deformed
sns.set(context='paper', style='white')
threeddefmatch = deftomatch.query('deformation=="Structural"')
deffinalmatch = datamatch_merge(
    threeddefmatch.query('type=="extrusion"'),
    threeddefmatch.query('type=="true-3D"'),
    "threshold",
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns="type")

# apply pe to all rows of dataframe matched, with threshold3d as the correct value and threshold as the estimated value
deffinalmatch["pe"] = deffinalmatch.apply(
    lambda row: pe_noabs(row["threshold3d"], row["threshold"], doabs=False), axis=1
)
plt.figure()
sns.boxplot(
    data=deffinalmatch.query('fiber_diam in [3,13]'),
    x="nerve_label",
    y="pe",
    hue="fiber_diam",
    # errorbar="se",
    # split=True,
    palette=sns.color_palette("RdPu", 2),
    whis=100,
    # width=1,
    # inner="point",
    # inner_kws=dict(box_width=2, whis_width=1)
)
# sns.stripplot(
#     data=deffinalmatch.query('fiber_diam in [3,13]'),
#     x="nerve_label",
#     y="pe",
#     hue="fiber_diam",
#     # errorbar="se",
#     # split=True,
#     palette=sns.color_palette("RdPu", 2),
#     # whis=100,
#     s=2,
#     edgecolor='k',
#     linewidth=0.3,
#     dodge=True
# )
# plt.ylim(-100, 100)
# plt.title('Threshold Percent Error by sample and fiber diameter')
plt.legend(title="D (μm)", ncols=3)
plt.xlabel("")
plt.ylabel("Relative Percent Difference (%)")
plt.gcf().set_size_inches([6, 5])
plt.axhline(0, color='black', ls='--', alpha=0.5)

# calculate min, max, and mean percent error for each fiber diameter
pe_means = deffinalmatch.groupby(["fiber_diam"]).agg(np.mean)
pe_medians = deffinalmatch.groupby(["fiber_diam"]).agg(np.median)
pe_mins = deffinalmatch.groupby(["fiber_diam"]).agg(np.min)
pe_maxs = deffinalmatch.groupby(["fiber_diam"]).agg(np.max)
print("Percent Error by Fiber Diameter")
print("Mean: ", pe_means["pe"])
print("Median: ", pe_medians["pe"])
print("Min: ", pe_mins["pe"])
print("Max: ", pe_maxs["pe"])
# now do the same but for absolute error new column 'ae' is the absolute value of column 'pe'
deffinalmatch["ae"] = deffinalmatch["pe"].abs()
ae_means = deffinalmatch.groupby(["fiber_diam"]).agg(np.mean)
ae_medians = deffinalmatch.groupby(["fiber_diam"]).agg(np.median)
ae_mins = deffinalmatch.groupby(["fiber_diam"]).agg(np.min)
ae_maxs = deffinalmatch.groupby(["fiber_diam"]).agg(np.max)
print("Absolute Error by Fiber Diameter")
print("Mean: ", ae_means["ae"])
print("Median: ", ae_medians["pe"])
print("Min: ", ae_mins["ae"])
print("Max: ", ae_maxs["ae"])
plt.gcf().set_size_inches(3, 2.5)
plt.figure()

sns.set(context='paper', style='white')
threeddefmatch = deftomatch.query('deformation=="Structural"')
deffinalmatch = datamatch_merge(
    threeddefmatch.query('type=="extrusion"'),
    threeddefmatch.query('type=="true-3D"'),
    "threshold",
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns="type")

# apply pe to all rows of dataframe matched, with threshold3d as the correct value and threshold as the estimated value
deffinalmatch["pe"] = deffinalmatch.apply(
    lambda row: pe_noabs(row["threshold3d"], row["threshold"], doabs=False), axis=1
)
plt.figure()
sns.barplot(
    data=deffinalmatch.query('fiber_diam in [3,13]'),
    x="nerve_label",
    y="pe",
    hue="fiber_diam",
    # errorbar="se",
    # split=True,
    palette=sns.color_palette("RdPu", 2),
    estimator='median',
    # errorbar=('sd',1)
    # whis=100,
    # width=1,
    # inner="point",
    # inner_kws=dict(box_width=2, whis_width=1)
)
# sns.stripplot(
#     data=deffinalmatch.query('fiber_diam in [3,13]'),
#     x="nerve_label",
#     y="pe",
#     hue="fiber_diam",
#     # errorbar="se",
#     # split=True,
#     palette=sns.color_palette("RdPu", 2),
#     # whis=100,
#     s=2,
#     edgecolor='k',
#     linewidth=0.3,
#     dodge=True
# )
# plt.ylim(-100, 100)
# plt.title('Threshold Percent Error by sample and fiber diameter')
plt.legend(title="D (μm)", ncols=3)
plt.xlabel("")
plt.ylabel("Relative Percent Difference (%)")
plt.gcf().set_size_inches([6, 5])
plt.axhline(0, color='black', ls='--', alpha=0.5)
plt.gcf().set_size_inches(3, 2.5)

# %% dose-response example
sns.set(context='paper', style='white', font_scale=1)
plt.figure()
drthis = defdr.query(
    f"fiber_diam in [3] and contact in {cath_comparison} and nerve_label =='3R' and deformation=='Structural'"
)

# Query the rows with type '3D'
df_3d = drthis.query("modeltype == 'true-3D'")

# Query the rows with type '2D' and sample ending in '1'
df_2d = drthis.query("modeltype == 'extrusion'")

# Merge the 3D and 2D data, keeping track of original row indices
merged_df = pd.merge(
    df_3d,
    df_2d,
    on=["nerve_label", "master_fiber_index", "nsim", "deformation"],
    suffixes=("_3d", "_2d"),
    how="left",
)  # TODO remove this how

# Update the 'inner', 'outer', and 'fiber' columns in the original DataFrame
drthis.loc[df_3d.index, "inner"] = merged_df["inner_2d"].values
drthis.loc[df_3d.index, "outer"] = merged_df["outer_2d"].values
drthis.loc[df_3d.index, "fiber"] = merged_df["fiber_2d"].values

g = sns.scatterplot(
    data=drthis,
    y="percent_activated",
    x="threshold",
    hue="inner",
    palette="rainbow",
    linewidth=0,
    s=50,
    marker=r"$\backslash$",
    # legend=False
)

sns.lineplot(
    data=drthis,
    y="percent_activated",
    x="threshold",
    color="k",
    linewidth=2,
    estimator=None,
    style="modeltype",
    legend=False,
    zorder=1,
    alpha=1,
)


plt.xlabel("Threshold (mA)")
plt.ylim([0, 1])
# plt.gcf().set_size_inches(8,4)
plt.ylabel("Proportion Activated")
# create legend, circle = extrusion, X= true-3D
# create handles manually

legend_elements = [
    Line2D([0], [0], label="extrusion", color="k", linestyle="-", alpha=0.6),
    Line2D([0], [0], label="true-3D", color="k", linestyle="--", alpha=0.6),
]
legend_labels = ["extrusion", "true-3D"]
g.legend(handles=legend_elements, labels=legend_labels, loc="lower right")
plt.gcf().set_size_inches(2, 2)

plt.figure()
g = sns.swarmplot(
    data=drthis,
    y="threshold",
    x="modeltype",
    hue="inner",
    palette="rainbow",
    linewidth=0,
    s=2,
    # legend=False
)
plt.xlabel("")
plt.ylabel("Threshold (mA)")
plt.legend([], [], frameon=False)
plt.ylim(0, None)
thismatch = datamatch_merge(
    drthis.query('modeltype=="extrusion"'),
    drthis.query('modeltype=="true-3D"'),
    "threshold",
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns="modeltype")
plt.gcf().set_size_inches(2, 2)

plt.figure()
sns.scatterplot(data=thismatch, x="threshold", y="threshold3d", hue="inner", palette="rainbow", legend=False, s=10)
plt.xlim(0, 4)
plt.ylim(0, 4)
plt.gca().set_aspect("equal")
plt.plot([0, 5], [0, 5], "k--")
plt.ylabel("True-3D Threshold")
plt.xlabel("Extrusion threshold")
plt.gcf().set_size_inches(2, 2)
# %% activation order fiberdiam
newdefdr = defdr.copy()
newdefdr = newdefdr.query("'6R' in nerve_label and deformation!='ASCENT'").sort_values("modeltype")
newdefdr["deformation"] = pd.Categorical(newdefdr["deformation"], categories=["Undeformed", "Structural"])
sns.set(font_scale=1.5)
sns.set_style("whitegrid")
newdefdr["percent_activated3"] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newdefdr.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        f"fiber_diam == 3 and nerve_label == @row.nerve_label and modeltype == @row.modeltype and sim == @row.sim and master_fiber_index == @row.master_fiber_index and contact in {cath_comparison} and deformation==@row.deformation"
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert val is not np.nan
    newdefdr.loc[row.Index, "percent_activated3"] = val
# %% plot
sns.set(font_scale=1, style="white", context='paper')
plt.figure()
newdefdr['fiber_diam'] = newdefdr['fiber_diam'].astype(str)
g = sns.catplot(
    kind="swarm",
    # row='deformation',
    data=newdefdr.query(f"fiber_diam in ['3','13'] and contact in {cath_comparison} and deformation=='Structural'"),
    y="percent_activated",
    x="fiber_diam",
    units="nerve_label",
    col="modeltype",
    palette="plasma",
    hue="percent_activated3",
    # hue='inner',
    estimator=None,
    linewidth=0,
    margin_titles=True,
    s=10,
)
xlim = plt.xlim()
# g.map_dataframe(sns.lineplot,
#                 y="percent_activated",
#                 x="fiber_diam",
#                 # hue='inner',
#                 estimator=None,
#                 units = 'master_fiber_index',
#                 color='black',
#                 alpha=0.3
#                 )
newdefdr['fiber_diam'] = newdefdr['fiber_diam'].astype(int)

plt.subplots_adjust(top=0.87)
# plt.suptitle(stringdat, x=0.37)
g.set_titles(row_template="", col_template="{col_name}")
g.axes[0][0].set_xlabel("")
g.axes[0][0].set_ylabel("Proportion fibers activated")
g.set_xlabels("D (μm)")
g.legend.remove()
norm = plt.Normalize(0, 1)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])

# Remove the legend and add a colorbar
g.figure.colorbar(
    sm,
    ax=g.axes.ravel().tolist(),
    aspect=10,
    shrink=0.8,
    label="Proportion 3 μm\nfibers activated",
    pad=0.1,
).ax.yaxis.set_ticks_position("left")

for i, mod in enumerate(newdefdr.modeltype.sort_values().unique()):
    shortdat = newdefdr.query(
        f"modeltype==@mod and fiber_diam in [13] and contact in {cath_comparison} and deformation=='Structural'"
    )
    data2d = shortdat.sort_values("percent_activated").master_fiber_index
    data3d = shortdat.sort_values("percent_activated3").master_fiber_index
    rc = compute_reorder_cost(list(data2d), list(data3d))
    ax = g.axes[0][i]
    text = ax.get_title() + f" - AR: {round(rc,3)}"
    ax.set_title(text)
plt.gcf().set_size_inches([5, 2])
plt.xlim(xlim)

# %% organizatino compare type and def
sns.set(font_scale=1, style='whitegrid', context='paper')
threshdat = newdefdat.copy()

scores = []
for nerve in pd.unique(threshdat["nerve_label"]):
    for deformtype in ["Undeformed", "Structural"]:
        for contact in ["3D", "cathodic", "center", "anodic"]:
            shortdat = threshdat.query(f'nerve_label=="{nerve}" and deformation ==@deformtype and contact in @contact')
            datasmol = shortdat.query("nsim==0").sort_values("threshold").master_fiber_index
            databeeg = shortdat.query("nsim==5").sort_values("threshold").master_fiber_index
            rc = compute_reorder_cost(list(datasmol), list(databeeg))
            scores.append(
                {
                    "sample": nerve,
                    "scoresmolbeeg": rc,
                    "deformation": deformtype,
                    "slice": contact,
                }
            )
fig = plt.figure()
scoredat = pd.DataFrame(scores)
sns.stripplot(
    data=scoredat,
    y="scoresmolbeeg",
    x="slice",
    color="black",
    hue="deformation",
    dodge=True,
    linewidth=0.5,
    edgecolor="white",
    s=4,
)
sns.barplot(
    data=scoredat, y="scoresmolbeeg", x="slice", hue="deformation", palette="binary", errorbar=None, estimator='median'
)
plt.ylabel("Activation Reordering")
plt.ylim(0, 0.6)
# plt.xlim(-1,2)
# plt.gcf().set_si
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(handles=handles[2:], labels=labs[2:])
# barplot
plt.xticks(
    [0, 1, 2, 3],
    ["true-3D", "extrusion\n(cathodic)", "extrusion\n(center)", "extrusion\n(anodic)"],
)
plt.axhline(0.45, color='red', linestyle='--')
plt.text(1.4, 0.47, 'random', color='red')
plt.xlabel("model")
plt.gcf().set_size_inches([3, 2])
# %% Percent Error deformed
sns.set(font_scale=1, style="white", context='paper')
# apply pe to all rows of dataframe matched, with threshold3d as the correct value and threshold as the estimated value
deffinalmatch["pe"] = deffinalmatch.apply(lambda row: pe(row["threshold3d"], row["threshold"]), axis=1)
plt.figure()
sns.boxplot(
    data=deffinalmatch.query('fiber_diam in [3,13]'),
    x="nerve_label",
    y="pe",
    hue="fiber_diam",
    # errorbar="se",
    # split=True,
    palette=sns.color_palette("RdPu", 2),
    whis=100,
    # width=1,
    # inner="point",
    # inner_kws=dict(box_width=2, whis_width=1)
)
# sns.stripplot(
#     data=deffinalmatch.query('fiber_diam in [3,13]'),
#     x="nerve_label",
#     y="pe",
#     hue="fiber_diam",
#     # errorbar="se",
#     # split=True,
#     palette=sns.color_palette("RdPu", 2),
#     # whis=100,
#     # width=1,
#     # inner="point",
#     # inner_kws=dict(box_width=2, whis_width=1)
# )
# plt.title('Threshold Percent Error by sample and fiber diameter')
legend = plt.legend(title="D (μm)", ncols=3)
plt.xlabel("")
plt.axhline(0, color='black', ls='--', alpha=0.5)

plt.ylabel("Absolute Percent Difference (%)")
legend.set_title("D (μm)")
plt.gcf().set_size_inches(3, 2.5)

plt.figure()
sns.barplot(
    data=deffinalmatch.query('fiber_diam in [3,13]'),
    x="nerve_label",
    y="pe",
    hue="fiber_diam",
    # errorbar="se",
    palette=sns.color_palette("RdPu", 2),
    # whis=100,
    estimator='median',
    # inner="point",
)
# plt.title('Threshold Percent Error by sample and fiber diameter')
legend = plt.legend(title="D (μm)", ncols=3, loc='lower left')
plt.xlabel("")
plt.axhline(0, color='black', ls='--', alpha=0.5)

plt.ylabel("Absolute Percent Difference (%)")
legend.set_title("D (μm)")
plt.gcf().set_size_inches(3, 2.5)
# plt.ylim(None,100)
# %% CCC comparison maincomp
nsimdata = concats.query("fiber_diam in [3,13]")  # TODO replace all cath comparison with non
g = sns.relplot(
    data=nsimdata.query("deformation=='Structural' and contact in @cath_comparison"),
    kind="scatter",
    col='nerve_label',
    x="threshold",
    y="threshold3d",
    # hue='Sample',
    row='fiber_diam',
    color="white",
    # s=20,
    facet_kws={"sharex": 'row', "sharey": 'row', "margin_titles": True},
    edgecolor="black",
    linewidth=1,
    alpha=1,
)
# for each plot, calculate CCC, and add to title
for pos, data in g.facet_data():
    nerve_label = data['nerve_label'].iloc[0]
    fiber_diam = data['fiber_diam'].iloc[0]
    for nerve in pd.unique(data["nerve_label"]):
        shortdat = data.query(f'nerve_label=="{nerve}"')
        data2d = shortdat["threshold"]
        data3d = shortdat["threshold3d"]
        ccc = concordance_correlation_coefficient(data3d, data2d)
        ax = g.axes[pos[0]][pos[1]]
        text = ax.get_title() + f" - CCC: {round(ccc,3)}"
        ax.set_title(text)
# for each row, find the max ylim and xlim, then set all axes to that, and plot a 1:1 line
for row in g.axes:
    maxlim = []
    for ax in row:
        maxlim.append(ax.get_xlim()[1])
        maxlim.append(ax.get_ylim()[1])
    maxlim = [0, max(maxlim)]
    for ax in row:
        ax.set_xlim(maxlim)
        ax.set_ylim(maxlim)
        ax.plot([0, maxlim[1]], [0, maxlim[1]], "k--")
# plt.gcf().set_size_inches(6, 6)
g.set_xlabels("Extrusion Threshold (mA)")
g.set_ylabels("True-3D Threshold (mA)")


# %% MCT
addln = True
imdata = pd.read_csv("thresh_unmatched_sim121_MCT.csv")
imdata = addpwfd(imdata, "121", infile="plotconfig_MCT")
imdata["type"] = imdata["type"].replace({"2D": "extrusion", "3D": "true-3D"})
imdata["nerve_label"] = imdata["nerve_label"].replace({"2L_MCT": "2L", "3R_MCT": "3R", "5R_MCT": "5R", "6R_MCT": "6R"})
imdata["active_src_index"] = imdata["active_src_index"].astype(str)
if (
    addln
):  # TODO figure out how to get this comparison working  - looks like centers of the two are slightly off and thus causing a disagree in fascicles?
    lndata = newdefdat.query('contact in @cath_comparison and deformation=="Structural" and fiber_diam in [3,13]')
    contname = "LN"
    lndata["active_src_index"] = contname
    imdata = pd.concat([imdata, lndata])
    if False:  # temp block which removes all but livanova
        imdata = imdata[imdata["active_src_index"] == contname]

imdata["contact"] = imdata["sample"].astype(str).str[2].replace({"2": "cathodic"})

imdata["fiber_diam"] = imdata["fiber_diam"].astype(int)

imdata.reset_index(inplace=True)

# inners for the true-3D data need to match their extrusion counterpart
for row in imdata.query('type=="true-3D"').itertuples():
    inner = imdata.query(
        f'type=="extrusion" and fiber_diam=={row.fiber_diam} and master_fiber_index=={row.master_fiber_index} and active_src_index=="{row.active_src_index}" and nerve_label=="{row.nerve_label}"'
    )["inner"]
    assert len(inner) == 1
    imdata.loc[row.Index, "inner"] = inner.values[0]

im_matched = datamatch_merge(
    imdata.query('type=="extrusion"'),
    imdata.query('type=="true-3D"'),
    "threshold",
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns="type")
# %% MCT CCC
sns.set(context="paper", font_scale=1, style="whitegrid")

g = sns.relplot(
    data=im_matched,
    kind="scatter",
    col="active_src_index",
    x="threshold",
    y="threshold3d",
    hue="nerve_label",
    row="fiber_diam",
    color="white",
    # s=20,
    palette=defpal,
    facet_kws={"sharex": False, "sharey": "row", "margin_titles": True},
    edgecolor="black",
    linewidth=1,
    alpha=1,
    # hue='inner',
)
# now calculate again but separate by nerve label
elceeceecee = []
# get data for each facet and calculate concordance correlation coefficient
for data in g.facet_data():
    for nerve in pd.unique(data[1]["nerve_label"]):
        shortdat = data[1].query(f'nerve_label=="{nerve}"')
        data2d = shortdat["threshold"]
        data3d = shortdat["threshold3d"]
        ccc = concordance_correlation_coefficient(data3d, data2d)
        elceeceecee.append(
            dict(
                contact=shortdat.active_src_index.unique()[0],
                fiber_diam=shortdat.fiber_diam.unique()[0],
                value=ccc**2,
                sample=nerve,
            )
        )
elcdata = pd.DataFrame(elceeceecee)

plt.figure()
elcdata["contact"] = elcdata["contact"].replace({"1": "MCT", "2": "MCT", "3": "MCT", "0": "MCT"})
elcdata["contact"] = elcdata["contact"].replace({"MCT": "Multi-Contact Cuff", "LN": "LivaNova Cuff"})
g = sns.barplot(
    data=elcdata,
    hue="fiber_diam",
    y="value",
    x="contact",
    palette=sns.color_palette("RdPu", 2),
    # edgecolor='k',
    errorbar=("sd", 1),
)
sns.stripplot(
    data=elcdata,
    y="value",
    x="contact",
    hue='fiber_diam',
    dodge=True,
    edgecolor='black',
    linewidth=0.5,
    s=3,
    palette=sns.color_palette("RdPu", 2),
    jitter=0.25,
)

plt.ylabel("CCC")
plt.xlabel("")
plt.legend(title="D (μm)", loc="lower left", ncol=2)
plt.ylim(0, 1)
plt.gcf().set_size_inches(3, 1.5)
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title="D: μm", handles=handles[:2], labels=labs[:2], ncol=2)
# %% MCT activation order
# dose response
imdr = imdata.copy()
imdr["active_src_index"] = pd.Categorical(
    imdata["active_src_index"],
    categories=sorted(imdata["active_src_index"].unique()),
    ordered=True,
)

imdr = calculate_dose_response(
    imdr,
    "threshold",
    "percent_activated",
    grouping_columns=["type", "sample", "fiber_diam", "sim", "active_src_index"],
)
imdr.sort_values("type", inplace=True)

# now do
newim = imdr.copy()
newim = newim.sort_values("type")
sns.set(font_scale=1.5)
sns.set_style("whitegrid")
newim["percent_activated2d"] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newim.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newim.query(
        'type == "true-3D" and fiber_diam == @row.fiber_diam and master_fiber_index == @row.master_fiber_index and active_src_index==@row.active_src_index and nerve_label==@row.nerve_label'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert val is not np.nan
    newim.loc[row.Index, "percent_activated3d"] = val
newim["type"] = pd.Categorical(newim["type"], ordered=True, categories=["true-3D", "extrusion"])

# %% plot
sns.set(font_scale=1, context='paper', style="whitegrid")
scores = []
for cc in sorted(imdata["active_src_index"].unique()):
    shortdat = imdata.query(f'active_src_index=="{cc}"')
    for nerve_label in pd.unique(shortdat["nerve_label"]):
        for fiber_diam in shortdat.fiber_diam.unique():
            finallythedat = shortdat.query(f'nerve_label=="{nerve_label}" and fiber_diam==@fiber_diam')
            data2d = finallythedat.query('type=="extrusion"').sort_values("threshold").master_fiber_index
            data3d = finallythedat.query('type=="true-3D"').sort_values("threshold").master_fiber_index
            rc = compute_reorder_cost(list(data2d), list(data3d))
            scores.append(
                {
                    "active_src_index": cc,
                    "score2d3d": rc,
                    "nerve_label": nerve_label,
                    "fiber_diam": fiber_diam,
                }
            )
scoredat = pd.DataFrame(scores)
scoredatnew = scoredat.copy()
scoredatnew["active_src_index"] = scoredatnew["active_src_index"].replace(
    {"1": "MCT", "2": "MCT", "3": "MCT", "0": "MCT"}
)
scoredatnew["active_src_index"] = scoredatnew["active_src_index"].replace(
    {"MCT": "Multi-Contact Cuff", "LN": "LivaNova Cuff"}
)
plt.figure()
g = sns.barplot(
    data=scoredatnew,
    y="score2d3d",
    x="active_src_index",
    hue="fiber_diam",
    palette=sns.color_palette("RdPu", 2),
    # edgecolor='k',
    errorbar=('sd', 1),
)
sns.stripplot(
    data=scoredatnew,
    y="score2d3d",
    x="active_src_index",
    hue='fiber_diam',
    dodge=True,
    edgecolor='black',
    linewidth=0.5,
    s=3,
    palette=sns.color_palette("RdPu", 2),
    jitter=0.25,
)

plt.xlabel("")
plt.ylabel("Activation Reordering")
plt.legend(title="D (μm)", ncol=2)
plt.gcf().set_size_inches(3, 1.5)
plt.ylim(0, 0.6)
plt.axhline(0.45, color='red', ls='--')
plt.text(-0.3, 0.5, 'random', color='red')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title="D: μm", handles=handles[:2], labels=labs[:2], ncol=2, loc='upper right')
# %% dose response MCT
sns.set(style="white", context='paper', font_scale=1)
imdr = imdata.copy().rename(columns={"sample": "samplenum", "type": "modeltype"})
imdr = calculate_dose_response(
    imdr,
    "threshold",
    "percent_activated",
    grouping_columns=[
        "modeltype",
        "samplenum",
        "fiber_diam",
        "sim",
        "active_src_index",
        "nerve_label",
    ],
)
imdr.sort_values("modeltype", inplace=True)
imdrmatch = datamatch_merge(
    imdr.query('modeltype=="extrusion"'),
    imdr.query('modeltype=="true-3D"'),
    "percent_activated",
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns="modeltype")

alldrcomp = []

for nerve_label in ['3R']:
    # get colors for size
    nerve_label_to_samplenum = {
        "3R": 3721,
    }
    samplenum = nerve_label_to_samplenum[nerve_label]
    # for samplenum,prefix in zip([3721,],['3R']):
    criteria = {
        'partial_matches': True,
        'include_downstream': False,
        'indices': {'sample': [samplenum], 'model': None, 'sim': None},
    }

    q = Query(criteria)
    q.run()

    results = q.summary()

    sample_index = results['samples'][0]['index']

    fig, ax = plt.subplots(1, 1)
    item: Sample = q.get_object(Object.SAMPLE, [results['samples'][0]['index']])
    slide = item.slides[0]
    cmap = cm.get_cmap('cmr.gem', len(slide.inners())).reversed()
    # plot fascicle ECD

    ECD = [{'inner': i, "D": f.ecd()} for i, f in enumerate(slide.inners())]
    ECD = pd.DataFrame(ECD)
    ECD.sort_values('D', inplace=True)
    colors = [to_hex(cmap(i)) for i in range(len(slide.inners()))]
    ECD['color'] = None
    for i, row in enumerate(ECD.itertuples()):
        ECD.loc[row.Index, 'color'] = colors[i]
    ECD.sort_values('inner', inplace=True)

    thisplotdata = imdr.query(f'nerve_label=="{nerve_label}"')
    fig, axs = plt.subplots(
        1, len(pd.unique(thisplotdata.active_src_index)), sharey=True, sharex=True
    )  # change to 5 to get LN
    for acsrc, ax in zip(sorted(pd.unique(thisplotdata.active_src_index)), axs):
        # plt.figure()
        g = sns.scatterplot(
            data=thisplotdata.query("fiber_diam in [3] and active_src_index==@acsrc"),
            y="percent_activated",
            x="threshold",
            hue="inner",
            palette=ECD['color'].values,
            hue_order=range(len(slide.fascicles)),
            linewidth=0,
            ax=ax,
            marker=r"$\backslash$",
            s=100,
            legend=False,
        )
        import matplotlib.patheffects as paffect

        sns.lineplot(
            data=thisplotdata.query("fiber_diam in [3] and active_src_index==@acsrc"),
            y="percent_activated",
            x="threshold",
            color="w",
            linewidth=1,
            estimator=None,
            units="modeltype",
            legend=False,
            zorder=1,
            ax=ax,
            style="modeltype",
            path_effects=[paffect.Stroke(linewidth=2.25, foreground='k'), paffect.Normal()],
        )

        ax.set_xlabel("Threshold (mA)")
        ax.set_ylim([0, 1])
        # create legend, circle = extrusion, X= true-3D
        # create handles manually
        from matplotlib.lines import Line2D

        legend_elements = [
            Line2D([0], [0], linestyle="-", color="k", label="extrusion", linewidth=2),
            Line2D([0], [0], linestyle="--", color="k", label="true-3D", linewidth=2),
        ]
        legend_labels = ["extrusion", "true-3D"]
        ax.set_title(f"Active contact: {acsrc}")
    g.legend(handles=legend_elements, labels=legend_labels, loc="lower right")
    axs[0].set_ylabel("Proportion fibers activated")
    plt.subplots_adjust(wspace=0)
    plt.gcf().set_size_inches(8, 1.5)

    # now do swarmplot
    plt.figure()
    g = sns.catplot(
        kind="swarm",
        data=imdata.query(f'fiber_diam==3 and nerve_label=="{nerve_label}"'),
        y="threshold",
        x="type",
        units="nerve_label",
        col="active_src_index",
        palette=ECD['color'].values,
        hue_order=range(len(slide.fascicles)),
        hue="inner",
        # hue='inner',
        estimator=None,
        linewidth=0,
        margin_titles=True,
        s=8,
        legend=False,
    )
    plt.subplots_adjust(top=0.87)
    # plt.suptitle(stringdat, x=0.37)
    g.set_titles(row_template="", col_template="Active contact: {col_name}")
    g.axes[0][0].set_xlabel("")
    g.axes[0][0].set_ylabel("Threshold (mA)")
    g.set_xlabels("")
    norm = plt.Normalize(0, 1)
    sm = plt.cm.ScalarMappable(cmap="rainbow", norm=norm)
    sm.set_array([])

    plt.gcf().set_size_inches(10, 3)
    plt.ylim(0, None)

    # plot slide

    lim = [-2000, 2000]
    max_thk = 1000
    r_cuff_in_pre_MCT = 1500
    criteria = {
        "partial_matches": True,
        "include_downstream": False,
        "indices": {"sample": [samplenum], "model": None, "sim": None},
    }

    q = Query(criteria)
    q.run()

    results = q.summary()

    sample_index = results["samples"][0]["index"]
    plt.figure()
    ax = plt.gca()

    item: Sample = q.get_object(Object.SAMPLE, [results["samples"][0]["index"]])
    slide = item.slides[0]
    cmap = cm.get_cmap("rainbow", len(slide.inners()))
    # plot a horizontal line across each sample
    # for each fascicle use orange (above line) for on target and blue (below line) for off target
    xcutoff = 0

    # add a circular band around the nerve with radius 1500 um from 0 to 270 degrees
    # theta = np.linspace(0, np.radians(325), 100)
    # r = 100+slide.nerve.ecd()/2
    R_in_MCT = slide.nerve.ecd() / 2 + 150
    thetamax = np.radians((1 * (r_cuff_in_pre_MCT / R_in_MCT)) * 360)
    theta = np.linspace(0, thetamax, 100)

    x = R_in_MCT * np.cos(theta)
    y = R_in_MCT * np.sin(theta)
    plt.plot(x, y, "k", linewidth=6, label="cuff")
    # add 6 evenly spaced stars from 30 to 300 degrees, add numbers 0-5 to the stars
    theta = np.linspace(np.radians(50), np.radians(325 - 50), 4)
    # add a circular band around the nerve with a certain arc length (e.g., 10 mm)

    contact_angs = []
    for i in range(4):
        ang_contactcenter_pre_MCT = 90
        ang_cuffseam_pre_MCT = 45
        ang_contactcenter_MCT = ang_contactcenter_pre_MCT * (r_cuff_in_pre_MCT / R_in_MCT)
        ang_cuffseam_MCT = ang_cuffseam_pre_MCT * (r_cuff_in_pre_MCT / R_in_MCT)
        contactpos = np.radians(ang_cuffseam_MCT + i * ang_contactcenter_MCT)
        contactspan = np.radians(15)
        contactheta = np.linspace(contactpos - contactspan, contactpos + contactspan, 100)
        contact_angs.append(contactpos)
        r = 110 + slide.nerve.ecd() / 2
        x = r * np.cos(contactheta)
        y = r * np.sin(contactheta)
        plt.plot(x, y, color="r", label="contacts" if i == 0 else "_", linewidth=4)
        txt = plt.text(
            x[int(len(x) / 2)],
            y[int(len(x) / 2)],
            str(i),
            fontsize=12,
            color="k",
            ha="center",
            va="center",
        )
        txt.set_path_effects([PathEffects.withStroke(linewidth=1.5, foreground="w")])

    slide.plot(
        fix_aspect_ratio=True,
        final=False,
        ax=ax,
        inner_index_labels=False,
        # scalebar=True,
        # scalebar_length=100,
        # scalebar_units='μm',
        fascicle_colors=ECD['color'].values,
        inner_format="w-",
        line_kws=dict(linewidth=1),
        colors_for_outers=True,
        inners_flag=False,
    )
    plt.xlabel("\u03bcm")
    plt.ylabel("\u03bcm")
    plt.ylim(lim)
    plt.xlim(lim)
    slide.add_scalebar(plt.gca(), 1, "mm")
    labels = ["on target", "off target"]
    plt.title(f"Nerve: {nerve_label}")
    # now calculate FASR
    imdata.reset_index(inplace=True, drop=True)
# %% calculate FASR complex
imdata.reset_index(inplace=True, drop=True)


def recruitment_cost(data, activated=1):
    """Calculate recruitment cost for each inner. :param activated: proportion
    of on-target fibers activated.

    Recruitment cost is defined as the ratio of number of stimulated
    off-target fibers to total number of off-target fibers. From
    https://iopscience.iop.org/article/10.1088/1741-2560/10/3/036010
    """
    for inner in pd.unique(data["inner"]):
        # get threshold for inner
        inner_data = data.query(f"inner == {inner}")["threshold"]
        assert len(data) > 200 & len(data) < 250
        # assert len(data['inner'])==len(set(data['inner']))
        # inner_thresh = np.amax(inner_data)
        # above line assumes 100% activation of on-target fibers, instead use activated
        inner_thresh = np.percentile(inner_data, activated * 100, method="higher")
        # get all off-target fiber thresholds
        off_thresh = data.query(f"inner != '{inner}'")["threshold"]
        # calculate recruitment cost
        cost = np.sum(off_thresh <= inner_thresh) / len(off_thresh)
        data.loc[data["inner"] == inner, "RC"] = cost
        # fascicle selectivity ratio is 1-RC
        data.loc[data["inner"] == inner, "FASR"] = 1 - cost
        fasr_dict = {
            "active_src_index": data["active_src_index"].iloc[0],
            "fiber_diam": data["fiber_diam"].iloc[0],
            "type": data["type"].iloc[0],
            "inner": inner,
            "RC": cost,
            "FASR": 1 - cost,
            "nerve_label": data["nerve_label"].iloc[0],
        }
        yield fasr_dict


imdatfasr = []
for contact_config in pd.unique(imdata["active_src_index"]):
    for fiber_diam in pd.unique(imdata["fiber_diam"]):
        for t in pd.unique(imdata["type"]):
            for nerve_label in pd.unique(imdata["nerve_label"]):
                imdatain = imdata.query(
                    f'active_src_index == "{contact_config}" and fiber_diam == {fiber_diam} and type == "{t}" and nerve_label=="{nerve_label}"'
                )
                imdatfasr.extend(recruitment_cost(imdatain, activated=0.10))

imdatfasr = pd.DataFrame(imdatfasr)

# %%
# rewrite the above with map dataframe
sns.set(style="whitegrid", context="paper", font_scale=1.5)
alldfcomp = pd.concat(alldrcomp)
alldfcomp["pe"] = alldfcomp.apply(lambda row: pe(row["threshold3d"], row["threshold"]), axis=1)
g = sns.FacetGrid(
    data=alldfcomp.query("fiber_diam==3"),
    palette=defpal,
    row="active_src_index",
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x="level",
    y="pe",
    marker="o",
    linewidth=1,
    markersize=5,
    alpha=0.6,
    hue="nerve_label",
    # estimator=None,
    palette=defpal,
)
g.map_dataframe(
    sns.lineplot,
    x="level",
    y="pe",
    marker="s",
    linewidth=2,
    markersize=6,
    estimator="median",
    errorbar=None,
    color="black",
)
plt.gcf().set_size_inches(3, 6)
plt.legend(ncol=2)
g.set_titles(col_template="{col_name}", row_template="")
rownames(g, row_template="Percent Error\n{row_name}")
for ax in g.axes.ravel():
    plt.sca(ax)
    plt.xlabel("Active Contact")
plt.ylim(None, 80)
# plt.suptitle('Black Square=median')
print(alldfcomp.groupby(['fiber_diam', 'active_src_index', 'level']).median())

# %% dose-response MCT final
sns.set_style("white")
imdr = imdata.copy().rename(columns={"sample": "samplenum", "type": "modeltype"})
imdr = calculate_dose_response(
    imdr,
    "threshold",
    "percent_activated",
    grouping_columns=[
        "modeltype",
        "samplenum",
        "fiber_diam",
        "sim",
        "active_src_index",
        "nerve_label",
    ],
)
imdr.sort_values("modeltype", inplace=True)
imdrmatch = datamatch_merge(
    imdr.query('modeltype=="extrusion"'),
    imdr.query('modeltype=="true-3D"'),
    "percent_activated",
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns="modeltype")

alldrcomp = []

for nerve_label in pd.unique(imdrmatch["nerve_label"]):
    thisplotdata = imdr.query(f'nerve_label=="{nerve_label}"')
    fig, axs = plt.subplots(
        1, len(pd.unique(thisplotdata.active_src_index)), sharey=True, sharex=True
    )  # change to 5 to get LN
    axs[0].set_ylabel("Proportion fibers activated\nNerve: " + nerve_label)
    plt.subplots_adjust(wspace=0)
    levels = {"onset": 10, "half": 50, "saturation": 90}
    grouped = thisplotdata.groupby(["nerve_label", "fiber_diam", "modeltype", "active_src_index"])
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
        id_vars=["active_src_index", "fiber_diam", "modeltype"],
        value_vars=["onset", "half", "saturation"],
        var_name="level",
        value_name="threshold",
    )

    # set up facetgrid with nsim as row and level as columns
    compiled_data.reset_index(inplace=True)
    # set fiber_diam to category
    compiled_data.modeltype = compiled_data.modeltype.astype("category")
    # add a units column with unique number for each combination of fiber_diam and level
    compiled_data["units"] = compiled_data.groupby(["fiber_diam", "level", "active_src_index"]).ngroup()
    compiled_data["fiber_diam"] = compiled_data["fiber_diam"].astype(int)
    compiled_data.dropna(inplace=True)

    compmatch = datamatch_merge(
        compiled_data.query('modeltype=="extrusion"'),
        compiled_data.query('modeltype=="true-3D"'),
        "threshold",
        merge_cols=["level", "active_src_index", "fiber_diam"],
    ).drop(columns="modeltype")

    compmatch["nerve_label"] = nerve_label

    alldrcomp.append(compmatch)

alldfcomp = pd.concat(alldrcomp)
alldfcomp["active_src_index"].replace({str(i): "MultiContact" for i in range(4)}, inplace=True)
sns.set(style="white", font_scale=1.5)
alldfcomp.rename(columns={"active_src_index": "cuff"}, inplace=True)

# rewrite the above with map dataframe
sns.set(style="white", context="paper", font_scale=1)
alldfcomp = pd.concat(alldrcomp)
alldfcomp["pe"] = alldfcomp.apply(lambda row: pe(row["threshold3d"], row["threshold"]), axis=1)
alldfcomp["active_src_index"] = alldfcomp["active_src_index"].replace({"1": "MCT", "2": "MCT", "3": "MCT", "0": "MCT"})
alldfcomp["active_src_index"] = alldfcomp["active_src_index"].replace(
    {"MCT": "Multi-Contact Cuff", "LN": "LivaNova Cuff"}
)
g = sns.FacetGrid(
    data=alldfcomp.query("fiber_diam==3"),
    palette=defpal,
    margin_titles=True,
)
g.map_dataframe(
    sns.pointplot,
    x="level",
    y="pe",
    marker="o",
    linewidth=1,
    markersize=5,
    alpha=0.6,
    hue="active_src_index",
    # estimator=None,
    palette='binary',
    # err_style='bars'
    dodge=True,
    estimator='median',
    errorbar=('pi', 25),
)

plt.gcf().set_size_inches(2, 2)
plt.legend()
g.set_titles(col_template="{col_name}", row_template="")
plt.xlabel('% active')
plt.xticks([0, 1, 2], ['10', '50', '90'])
plt.ylabel('Absolute Percent Difference (%)')
plt.ylim(None, 30)
# plt.suptitle('Black Square=median')
# %% dose-response MCT final
# rewrite the above with map dataframe
sns.set(style="white", context="paper", font_scale=1)
alldfcomp = pd.concat(alldrcomp)
alldfcomp["pe"] = alldfcomp.apply(lambda row: pe(row["threshold3d"], row["threshold"]), axis=1)
alldfcomp["active_src_index"] = alldfcomp["active_src_index"].replace({"1": "MCT", "2": "MCT", "3": "MCT", "0": "MCT"})
alldfcomp["active_src_index"] = alldfcomp["active_src_index"].replace(
    {"MCT": "Multi-Contact Cuff", "LN": "LivaNova Cuff"}
)
alldfcomp.fiber_diam = alldfcomp.fiber_diam.astype(str)
g = sns.FacetGrid(data=alldfcomp, palette=defpal, margin_titles=True, col='fiber_diam')
g.map_dataframe(
    sns.pointplot,
    x="level",
    y="pe",
    marker="s",
    linewidth=1,
    markersize=5,
    alpha=1,
    hue="active_src_index",
    # estimator=None,
    palette='binary',
    # err_style='bars'
    dodge=True,
    estimator='median',
    markeredgecolor='black',
)
g.map_dataframe(
    sns.swarmplot,
    x="level",
    y="pe",
    s=3,
    alpha=1,
    hue="active_src_index",
    # estimator=None,
    palette='binary',
    # err_style='bars'
    dodge=True,
    # zorder=-1,
    linewidth=0,
)

plt.gcf().set_size_inches(4, 2)
plt.legend()
g.set_titles(col_template="D: {col_name} μm", row_template="")
g.set_xlabels('% active')
for ax in g.axes.ravel():
    ax.tick_params(pad=0)
plt.xticks([0, 1, 2], ['10', '50', '90'])
plt.ylabel('Absolute Percent Difference (%)')
plt.ylim(None, 30)
g.set_ylabels('Absolute Percent Difference (%)')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title="", handles=handles[:2], labels=labs[:2], frameon=False)
# plt.suptitle('Black Square=median')
print(alldfcomp.groupby(['fiber_diam', 'active_src_index', 'level']).median())
# %% dose-response MCT final noabs
# rewrite the above with map dataframe
sns.set(style="white", context="paper", font_scale=1)
alldfcomp = pd.concat(alldrcomp)
alldfcomp["pe"] = alldfcomp.apply(lambda row: pe_noabs(row["threshold3d"], row["threshold"]), axis=1)
alldfcomp["active_src_index"] = alldfcomp["active_src_index"].replace({"1": "MCT", "2": "MCT", "3": "MCT", "0": "MCT"})
alldfcomp["active_src_index"] = alldfcomp["active_src_index"].replace(
    {"MCT": "Multi-Contact Cuff", "LN": "LivaNova Cuff"}
)
alldfcomp.fiber_diam = alldfcomp.fiber_diam.astype(str)
g = sns.FacetGrid(data=alldfcomp, palette=defpal, margin_titles=True, col='fiber_diam')
g.map_dataframe(
    sns.pointplot,
    x="level",
    y="pe",
    marker="s",
    linewidth=1,
    markersize=5,
    # alpha=0.6,
    hue="active_src_index",
    # estimator=None,
    palette='binary',
    # err_style='bars'
    dodge=True,
    estimator='median',
    markeredgecolor='black',
)
g.map_dataframe(
    sns.swarmplot,
    x="level",
    y="pe",
    marker="o",
    linewidth=0,
    # markersize=5,
    # alpha=0.6,
    s=3,
    hue="active_src_index",
    # estimator=None,
    palette='binary',
    # err_style='bars'
    # dodge=True,
    # estimator=None,
    # units = 'nerve_label'
    dodge=True,
)
for ax in g.axes.ravel():
    ax.axhline(0, color='k', ls='--')
for ax in g.axes.ravel():
    ax.tick_params(pad=0)
g.set_ylabels('Relative Percent Difference (%)')
plt.gcf().set_size_inches(4, 2)
# plt.legend()
g.set_titles(col_template="D: {col_name} μm", row_template="")
g.set_xlabels('% active')
plt.xticks([0, 1, 2], ['10', '50', '90'])
# plt.ylabel('Absolute Percent Difference (%)')
plt.ylim(None, 30)
# plt.suptitle('Black Square=median')
# %% FASR by location
imdata.reset_index(inplace=True, drop=True)


def recruitment_cost(data, activated=1):
    """Calculate recruitment cost for each inner. :param activated: proportion
    of on-target fibers activated.

    Recruitment cost is defined as the ratio of number of stimulated
    off-target fibers to total number of off-target fibers. From
    https://iopscience.iop.org/article/10.1088/1741-2560/10/3/036010
    """
    for inner in pd.unique(data["inner"]):
        # get threshold for inner
        inner_data = data.query(f"inner == {inner}")["threshold"]
        assert len(data) > 200 & len(data) < 250
        # assert len(data['inner'])==len(set(data['inner']))
        # inner_thresh = np.amax(inner_data)
        # above line assumes 100% activation of on-target fibers, instead use activated
        inner_thresh = np.percentile(inner_data, activated * 100, method="higher")
        # get all off-target fiber thresholds
        off_thresh = data.query(f"inner != '{inner}'")["threshold"]
        # calculate recruitment cost
        cost = np.sum(off_thresh <= inner_thresh) / len(off_thresh)
        data.loc[data["inner"] == inner, "RC"] = cost
        # fascicle selectivity ratio is 1-RC
        data.loc[data["inner"] == inner, "FASR"] = 1 - cost
        fasr_dict = {
            "active_src_index": data["active_src_index"].iloc[0],
            "fiber_diam": data["fiber_diam"].iloc[0],
            "type": data["type"].iloc[0],
            "inner": inner,
            "RC": cost,
            "FASR": 1 - cost,
            "nerve_label": data["nerve_label"].iloc[0],
        }
        yield fasr_dict


imdatfasr = []
for contact_config in pd.unique(imdata["active_src_index"]):
    for fiber_diam in pd.unique(imdata["fiber_diam"]):
        for t in pd.unique(imdata["type"]):
            for nerve_label in pd.unique(imdata["nerve_label"]):
                imdatain = imdata.query(
                    f'active_src_index == "{contact_config}" and fiber_diam == {fiber_diam} and type == "{t}" and nerve_label=="{nerve_label}"'
                )
                imdatfasr.extend(recruitment_cost(imdatain, activated=0.10))

imdatfasr = pd.DataFrame(imdatfasr)
# %% MCT selectivity
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(font_scale=1, style="white", context='paper')

lim = [-2000, 2000]
max_thk = 1000
analysisdiam = 3

# sns.set(font_scale=1,style='white')
for nerve_label, samplenum, r_cuff_in_pre_MCT in zip(
    ["2L", "3R", "5R", "6R"], [25212, 37212, 57212, 67212], [1000, 1500, 1000, 1000]
):
    fig, axs = plt.subplots(1, 2, width_ratios=(1, 1))
    plt.gcf().set_size_inches(10, 4)
    criteria = {
        "partial_matches": True,
        "include_downstream": False,
        "indices": {"sample": [samplenum], "model": None, "sim": None},
    }

    q = Query(criteria)
    q.run()

    results = q.summary()

    sample_index = results["samples"][0]["index"]
    ax = axs[0]
    plt.sca(ax)

    item: Sample = q.get_object(Object.SAMPLE, [results["samples"][0]["index"]])
    slide = item.slides[0]
    cmap = cm.get_cmap("rainbow", len(slide.inners()))
    # plot a horizontal line across each sample
    # for each fascicle use orange (above line) for on target and blue (below line) for off target
    xcutoff = 0

    # add a circular band around the nerve with radius 1500 um from 0 to 270 degrees
    # theta = np.linspace(0, np.radians(325), 100)
    # r = 100+slide.nerve.ecd()/2
    R_in_MCT = slide.nerve.ecd() / 2 + 150
    thetamax = np.radians((1 * (r_cuff_in_pre_MCT / R_in_MCT)) * 360)
    theta = np.linspace(0, thetamax, 100)

    x = R_in_MCT * np.cos(theta)
    y = R_in_MCT * np.sin(theta)
    plt.plot(x, y, "k", linewidth=6, label="cuff")
    # add 6 evenly spaced stars from 30 to 300 degrees, add numbers 0-5 to the stars
    theta = np.linspace(np.radians(50), np.radians(325 - 50), 4)
    # add a circular band around the nerve with a certain arc length (e.g., 10 mm)

    contact_angs = []
    for i in range(4):
        ang_contactcenter_pre_MCT = 90
        ang_cuffseam_pre_MCT = 45
        ang_contactcenter_MCT = ang_contactcenter_pre_MCT * (r_cuff_in_pre_MCT / R_in_MCT)
        ang_cuffseam_MCT = ang_cuffseam_pre_MCT * (r_cuff_in_pre_MCT / R_in_MCT)
        contactpos = np.radians(ang_cuffseam_MCT + i * ang_contactcenter_MCT)
        contactspan = np.radians(15)
        contactheta = np.linspace(contactpos - contactspan, contactpos + contactspan, 100)
        contact_angs.append(contactpos)
        r = 110 + slide.nerve.ecd() / 2
        x = r * np.cos(contactheta)
        y = r * np.sin(contactheta)
        plt.plot(x, y, color="r", label="contacts" if i == 0 else "_", linewidth=4)
        txt = plt.text(
            x[int(len(x) / 2)],
            y[int(len(x) / 2)],
            str(i),
            fontsize=12,
            color="k",
            ha="center",
            va="center",
        )
        txt.set_path_effects([PathEffects.withStroke(linewidth=1.5, foreground="w")])
    # add legend for on and off target
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="on target",
            markerfacecolor="turquoise",
            markersize=6,
            linewidth=0,
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="off target",
            markerfacecolor="red",
            markersize=6,
            linewidth=0,
        ),
    ]

    angle_radians = np.mean(contact_angs)

    # Calculate the slope of the line
    m = math.tan(angle_radians)

    ontarget = []
    offtarget = []
    targetcolors = []
    for i, fascicle in enumerate(slide.fascicles):
        centroid_x, centroid_y = fascicle.inners[0].centroid()

        # Check if the centroid is above the line y = m*x passing through the origin
        if centroid_y > m * centroid_x:
            targetcolors.append("turquoise")
            ontarget.append(i)
        else:
            targetcolors.append("red")
            offtarget.append(i)
    slide.plot(
        fix_aspect_ratio=True,
        final=False,
        ax=ax,
        inner_index_labels=False,
        # scalebar=True,
        # scalebar_length=100,
        # scalebar_units='μm',
        fascicle_colors=targetcolors,
        inner_format="w-",
        line_kws=dict(linewidth=1),
        colors_for_outers=True,
        inners_flag=False,
    )
    plt.xlabel("\u03bcm")
    plt.ylabel("\u03bcm")
    plt.ylim(lim)
    plt.xlim(lim)
    slide.add_scalebar(plt.gca(), 1, "mm")
    labels = ["on target", "off target"]
    plt.legend(handles=handles, labels=labels, loc="lower left", ncol=1, bbox_to_anchor=[0, -0.2])
    plt.title(f"Nerve: {nerve_label}")
    x_values = np.linspace(
        np.amin(slide.nerve.points[:, 0]) * 2, 2 * np.amax(slide.nerve.points[:, 0]), 400
    )  # Adjust range as needed
    y_values = m * x_values
    plt.plot(x_values, y_values, linestyle="--", color="k")
    # now calculate FASR
    imdata.reset_index(inplace=True, drop=True)

    def recruitment_cost(data, activated=1, targetcol="inner"):
        """Calculate recruitment cost for each inner. :param activated:
        proportion of on-target fibers activated.

        Recruitment cost is defined as the ratio of number of stimulated
        off-target fibers to total number of off-target fibers. From
        https://iopscience.iop.org/article/10.1088/1741-2560/10/3/036010
        """
        for inner in pd.unique(data[targetcol]):
            # get threshold for inner
            inner_data = data.query(f"{targetcol} == '{inner}'")["threshold"]
            assert len(data) > 200 & len(data) < 250
            # assert len(data['inner'])==len(set(data['inner']))
            # inner_thresh = np.amax(inner_data)
            # above line assumes 100% activation of on-target fibers, instead use activated
            inner_thresh = np.percentile(inner_data, activated * 100, method="higher")
            # get all off-target fiber thresholds
            off_thresh = data.query(f"{targetcol} != '{inner}'")["threshold"]
            # calculate recruitment cost
            cost = np.sum(off_thresh <= inner_thresh) / len(off_thresh)
            data.loc[data[targetcol] == inner, "RC"] = cost
            # fascicle selectivity ratio is 1-RC
            data.loc[data[targetcol] == inner, "FASR"] = 1 - cost
            fasr_dict = {
                "active_src_index": data["active_src_index"].iloc[0],
                "fiber_diam": data["fiber_diam"].iloc[0],
                "type": data["type"].iloc[0],
                targetcol: inner,
                "RC": cost,
                "FASR": 1 - cost,
                "nerve_label": data["nerve_label"].iloc[0],
                "percent_ontarget": activated,
            }
            yield fasr_dict

    imdatfasr = []
    for contact_config in pd.unique(imdata["active_src_index"]):
        for fiber_diam in pd.unique(imdata["fiber_diam"]):
            for t in pd.unique(imdata["type"]):
                for percent_ontarget in [0.1, 0.5, 0.9]:
                    imdatain = imdata.query(
                        f'active_src_index == "{contact_config}" and fiber_diam == {fiber_diam} and type == "{t}" and nerve_label=="{nerve_label}"'
                    )
                    imdatain["topo"] = imdatain.apply(
                        lambda row: ("on_target" if row["inner"] in ontarget else "off_target"),
                        axis=1,
                    )
                    imdatain["percent_ontarget"] = percent_ontarget
                    imdatfasr.extend(recruitment_cost(imdatain, activated=percent_ontarget, targetcol="topo"))
    imdatfasr = pd.DataFrame(imdatfasr)
    # now new plot, copy as above where extrusion and true-3D are plotted on the same graph with different markers, and the pointplot is added
    ax = axs[1]
    plt.sca(ax)
    imdatplot = imdatfasr.query(f'fiber_diam in [{analysisdiam}] and topo=="on_target"')
    imdatplot['RC'] *= 100

    g = sns.stripplot(
        data=imdatplot.query('type=="true-3D"'),
        x="active_src_index",
        y="RC",
        hue="percent_ontarget",
        palette="plasma",
        legend=False,
        marker="o",
        dodge=True,
        edgecolor="black",
        linewidth=1,
        s=5,
    )
    g = sns.stripplot(
        data=imdatplot.query('type=="extrusion"'),
        x="active_src_index",
        y="RC",
        hue="percent_ontarget",
        # palette="plasma",
        legend=False,
        marker="x",
        dodge=True,
        # edgecolor="black",
        linewidth=1,
        color='black',
        s=5,
    )
    sns.pointplot(
        data=imdatplot,
        x="active_src_index",
        y="RC",
        hue="percent_ontarget",
        palette="plasma",
        legend=False,
        marker=None,
        err_kws={"linewidth": 2},
        dodge=0.52,
        linestyle="none",
        errorbar=("pi", 100),
    )
    plt.ylabel("Percent off-target activated")
    # vertical line dashed in between each x value
    for i in range(5 if addln else 4):
        plt.axvline(i + 0.5, color="black", ls="--")
    # add legend for the two types. Create handles manually (gray marker with black outline)
    # also add line elements for each of the 4 percent_ontarget values from plasma colormap
    plasma = cm.get_cmap("plasma", 3)

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="true-3D",
            markerfacecolor="gray",
            # markersize=10,
            markeredgewidth=1,
            markeredgecolor="black",
        ),
        Line2D(
            [0],
            [0],
            marker="x",
            color="w",
            label="extrusion",
            markerfacecolor="black",
            # markersize=10,
            markeredgewidth=1,
            markeredgecolor="black",
        ),
        Line2D(
            [0],
            [0],
            linestyle="-",
            color=plasma(0),
            label="% on target fibers active",
            linewidth=0,
        ),
        Line2D([0], [0], linestyle="-", color=plasma(0), label="10%", linewidth=2),
        Line2D([0], [0], linestyle="-", color=plasma(1), label="50%", linewidth=2),
        Line2D([0], [0], linestyle="-", color=plasma(2), label="90%", linewidth=2),
    ]
    legend_labels = [
        "true-3D",
        "extrusion",
        "\ntarget fibers active",
        "10% (onset)",
        "50% (half)",
        "90% (saturation)",
    ]
    plt.legend(handles=legend_elements, labels=legend_labels, loc=(1.05, 0.2))
    plt.ylabel("off-target activated (%)")
    plt.ylim(-10, 110)
    plt.xlabel("Active contact")
    if addln:
        plt.xticks(range(5), list(range(4)) + ["LN"])
    else:
        plt.xticks(range(4), list(range(4)))
    plt.title(f"Nerve: {nerve_label} - D: {analysisdiam} μm")
    axs[0].axis('off')
    axs[0].set_title('')
    axs[1].set_title(f'D: {analysisdiam} μm')
    plt.gcf().set_size_inches(6, 2.5)
    # savef(f'MCT_{nerve_label} diam {analysisdiam}')
    # sys.exit()
    plt.figure()
# %% MCT selectivity just the numbers
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(font_scale=1, style="white", context='paper')
import matplotlib.patheffects as PathEffects

lim = [-2000, 2000]
max_thk = 1000
allselectdata = []
for analysisdiam in imdata.fiber_diam.unique():
    print(analysisdiam)
    # sns.set(font_scale=1,style='white')
    for nerve_label, samplenum, r_cuff_in_pre_MCT in zip(
        ["2L", "3R", "5R", "6R"], [25212, 37212, 57212, 67212], [1000, 1500, 1000, 1000]
    ):
        imdata.reset_index(inplace=True, drop=True)

        def recruitment_cost(data, activated=1, targetcol="inner"):
            """Calculate recruitment cost for each inner. :param activated:
            proportion of on-target fibers activated.

            Recruitment cost is defined as the ratio of number of stimulated
            off-target fibers to total number of off-target fibers. From
            https://iopscience.iop.org/article/10.1088/1741-2560/10/3/036010
            """
            for inner in pd.unique(data[targetcol]):
                # get threshold for inner
                inner_data = data.query(f"{targetcol} == '{inner}'")["threshold"]
                assert len(data) > 200 & len(data) < 250
                # assert len(data['inner'])==len(set(data['inner']))
                # inner_thresh = np.amax(inner_data)
                # above line assumes 100% activation of on-target fibers, instead use activated
                inner_thresh = np.percentile(inner_data, activated * 100, method="higher")
                # get all off-target fiber thresholds
                off_thresh = data.query(f"{targetcol} != '{inner}'")["threshold"]
                # calculate recruitment cost
                cost = np.sum(off_thresh <= inner_thresh) / len(off_thresh)
                data.loc[data[targetcol] == inner, "RC"] = cost
                # fascicle selectivity ratio is 1-RC
                data.loc[data[targetcol] == inner, "FASR"] = 1 - cost
                fasr_dict = {
                    "active_src_index": data["active_src_index"].iloc[0],
                    "fiber_diam": data["fiber_diam"].iloc[0],
                    "type": data["type"].iloc[0],
                    targetcol: inner,
                    "RC": cost,
                    "FASR": 1 - cost,
                    "nerve_label": data["nerve_label"].iloc[0],
                    "percent_ontarget": activated,
                }
                yield fasr_dict

        imdatfasr = []
        for contact_config in pd.unique(imdata["active_src_index"]):
            for fiber_diam in pd.unique(imdata["fiber_diam"]):
                for t in pd.unique(imdata["type"]):
                    for percent_ontarget in [0.1, 0.5, 0.9]:
                        imdatain = imdata.query(
                            f'active_src_index == "{contact_config}" and fiber_diam == {fiber_diam} and type == "{t}" and nerve_label=="{nerve_label}"'
                        )
                        imdatain["topo"] = imdatain.apply(
                            lambda row: ("on_target" if row["inner"] in ontarget else "off_target"),
                            axis=1,
                        )
                        imdatain["percent_ontarget"] = percent_ontarget
                        imdatfasr.extend(recruitment_cost(imdatain, activated=percent_ontarget, targetcol="topo"))
        imdatfasr = pd.DataFrame(imdatfasr)
        # now new plot, copy as above where extrusion and true-3D are plotted on the same graph with different markers, and the pointplot is added
        ax = axs[1]
        plt.sca(ax)
        imdatplot = imdatfasr.query(f'fiber_diam in [{analysisdiam}] and topo=="on_target"')
        imdatplot['RC'] *= 100
        imdatplot['fiber_diam'] = analysisdiam
        imdatplot['nerve_label'] = nerve_label
        allselectdata.append(imdatplot)
fiinalsel = datamatch_merge(
    pd.concat(allselectdata).query('type=="extrusion"'),
    pd.concat(allselectdata).query('type=="true-3D"'),
    "RC",
    merge_cols=["active_src_index", "fiber_diam", "nerve_label", "percent_ontarget"],
).drop(columns="type")
fiinalsel['resid'] = fiinalsel.RC3d - fiinalsel.RC
g = sns.FacetGrid(
    data=fiinalsel,
    row='fiber_diam',
)
g.map_dataframe(
    sns.stripplot,
    x='active_src_index',
    y='resid',
    hue='percent_ontarget',
    palette='plasma',
    legend=False,
    dodge=True,
    jitter=False,
    linewidth=1,
    alpha=0.6,
    zorder=-2,
)
g.map_dataframe(
    sns.pointplot,
    x='active_src_index',
    y='resid',
    hue='percent_ontarget',
    palette='plasma',
    legend=False,
    estimator='median',
    errorbar=None,
    markeredgewidth=1,
    markeredgecolor='k',
    marker='s',
)

plt.gcf().set_size_inches(2.5, 8)
plt.ylim(-100, 100)
g.set_titles(row_template='D: {row_name} μm')
g.set_ylabels('Off target activation (%) \n(true-3D minus extrusion)')
g.set_xlabels('Active Contact')
for ax in g.axes.ravel():
    for pos, col in zip(np.arange(0.5, 4, 1), ['gray', 'gray', 'gray', 'k']):
        ax.axvline(pos, ls='--', color=col)
# now absolute
fiinalsel['absresid'] = np.abs(fiinalsel.RC3d - fiinalsel.RC)
g = sns.FacetGrid(
    data=fiinalsel,
    row='fiber_diam',
)
g.map_dataframe(
    sns.barplot,
    x='active_src_index',
    y='absresid',
    hue='percent_ontarget',
    # palette='plasma',
    palette=['white'] * 3,
    edgecolor='k',
    legend=False,
    estimator='median',
    errorbar=None,
    # markeredgewidth=1,
    # markeredgecolor='k',
    # marker='s',
)
g.map_dataframe(
    sns.stripplot,
    x='active_src_index',
    y='absresid',
    hue='percent_ontarget',
    palette='plasma',
    legend=False,
    dodge=True,
    jitter=False,
    linewidth=1,
    alpha=0.6,
    # zorder=-2,
)
plt.gcf().set_size_inches(2.5, 4)
plt.ylim(0, 100)
g.set_titles(row_template='D: {row_name} μm')
g.set_ylabels('Off target activation (%) \n|true-3D minus extrusion|')
g.set_xlabels('Active Contact')
for ax in g.axes.ravel():
    for pos, col in zip(np.arange(0.5, 4 if addln else 3, 1), ['gray', 'gray', 'gray', 'k']):
        ax.axvline(pos, ls='--', color=col)
simipledat = fiinalsel.copy()
simipledat['active_src_index'] = simipledat['active_src_index'].replace(
    {'1': 'MCT', '2': 'MCT', '3': 'MCT', '0': 'MCT'}
)
print(simipledat.groupby(['fiber_diam', 'active_src_index'])['resid', 'absresid'].median())
print(simipledat.groupby(['fiber_diam', 'active_src_index'])['resid', 'absresid'].min())
print(simipledat.groupby(['fiber_diam', 'active_src_index'])['resid', 'absresid'].max())

# lastly, calculate the min off target activated for each nerve, fiber_diam, active_src_index, and percent_ontarget
selectivedat = simipledat.groupby(['nerve_label', 'fiber_diam', 'active_src_index', 'percent_ontarget'])['RC3d'].min()
# plot
g = sns.FacetGrid(
    data=selectivedat.reset_index(),
    row='fiber_diam',
)
g.map_dataframe(
    sns.barplot,
    x='active_src_index',
    y='RC3d',
    hue='percent_ontarget',
    palette='plasma',
    dodge=True,
    errorbar=None,
)
g.map_dataframe(
    sns.stripplot,
    x='active_src_index',
    color='k',
    y='RC3d',
    hue='percent_ontarget',
    dodge=True,
    marker='o',
    s=5,
    alpha=0.6,
    edgecolor='white',
    linewidth=1,
)
g.set_titles(row_template='D: {row_name} μm')
g.set_ylabels('Min. off target (%)')
g.set_xlabels('Cuff')
for ax in g.axes.ravel():
    ax.set_xticks([0, 1], ['LivaNova', 'MultiContact'])
plt.gcf().set_size_inches(2.5, 4)
# for ax in g.axes.ravel(): #uncomment for numbers
#     for i in ax.containers:
#         ax.bar_label(i,fmt='{:.2f}')
# %% threshold variances and coefficient of variation intrafascicle and inter
estimator, errorbar = "median", ('ci', 95)

sns.set(font_scale=1, style="white", context='paper')
vardat = repeated_deformation.query('contact=="cathodic" and deformation=="Structural"')
grouped = vardat.groupby(["contact", "fiber_diam", "type", "inner", "sample"])
analysis = grouped.agg({"threshold": [np.var, np.mean, variation]})
analysis.columns = ["_".join(col_name).rstrip("_") for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)


plt.figure()
g = sns.pointplot(
    data=analysis,
    y="threshold_variation",
    x="fiber_diam",
    hue="type",
    palette=pal2d3d,
    dodge=True,
    estimator=estimator,
    errorbar=errorbar,
)
plt.title("intrafascicle")
plt.ylabel("Threshold CoV")
plt.xlabel("Fiber Diameter (μm)")
# plt.yscale('log')
plt.gca().get_legend().set_title("")
plt.gcf().set_size_inches([6, 4])
plt.ylim(0, None)
plt.gcf().set_size_inches(3, 2)

# now do variance between fascicle mean thresholds
grouped = analysis.groupby(["contact", "fiber_diam", "type", "sample"])
analysis = grouped.agg({"threshold_mean": [np.var, variation, np.count_nonzero]})
analysis.columns = ["_".join(col_name).rstrip("_") for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)

plt.figure()
g = sns.pointplot(
    data=analysis,
    y="threshold_mean_variation",
    x="fiber_diam",
    hue="type",
    palette=pal2d3d,
    dodge=True,
    estimator=estimator,
    errorbar=errorbar,
)
plt.title("interfascicle")
plt.ylabel("Threshold CoV")
plt.xlabel("Fiber Diameter (μm)")
plt.gca().get_legend().set_title("")
plt.ylim(0, None)
plt.gcf().set_size_inches(3, 2)

# threshold variances for whole nerve
vardat = repeated_deformation.query('contact=="cathodic" and deformation=="Structural"')
grouped = vardat.groupby(["contact", "fiber_diam", "sim", "type", "sample"])
analysis = grouped.agg({"threshold": [np.var, np.mean, variation]})
analysis.columns = ["_".join(col_name).rstrip("_") for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)

plt.figure()
g = sns.pointplot(
    data=analysis,
    y="threshold_variation",
    x="fiber_diam",
    hue="type",
    palette=pal2d3d,
    dodge=True,
    estimator=estimator,
    errorbar=errorbar,
)
plt.title("intra-sample")
plt.ylabel("Threshold CoV")
plt.xlabel("Fiber Diameter (μm)")
plt.gca().get_legend().set_title("")
plt.ylim(0, None)
plt.gcf().set_size_inches(3, 2)

# %% 2d higher than 3d
sns.set(style="white", context="paper")
# do as facetgrid with col=type
drthis = defdr.query(
    f"fiber_diam in [3,13] and contact in {cath_comparison} and nerve_label =='3R' and deformation=='Structural'"
)

g = sns.catplot(
    kind='violin',
    data=defdr.query(f"fiber_diam in [3,13] and contact in {cath_comparison} and deformation=='Structural'"),
    x='nerve_label',
    hue="modeltype",
    palette=pal2d3d,
    y="threshold",
    split=True,
    col='fiber_diam',
    sharey=False,
    # facet_kws = {'sharex':'col'}
)
# plt.xlim(0, 4)
# plt.ylim(0, 1)
g.legend.set_title('')
plt.gcf().set_size_inches(8, 4)
plt.subplots_adjust(hspace=0.2)
g.axes[0][0].set_ylabel("Threshold (mA)")
g.set_titles(col_template="D: {col_name} μm")
plt.ylabel('')
g.set_xlabels('')
# %% residuals 2D3D
datacompare = matched_deformation.copy()
datacompare["residual"] = datacompare["threshold"] - datacompare["threshold3d"]
sns.set(style="white", context="paper")
g = sns.catplot(
    data=datacompare,
    kind='violin',
    y="residual",
    # split modeltype for violinplot
    # palette=pal2d3d,
    split=True,
    x=fiber_diam,
    col="nerve_label",
)


# %% threshold range
estimator, errorbar = "median", ('ci', 95)
from scipy.stats import variation


def rangecalc(data):
    return np.amax(data) - np.amin(data)


sns.set(font_scale=1, style="white", context='paper')
vardat = repeated_deformation.query('contact=="cathodic" and deformation=="Structural"')
grouped = vardat.groupby(["contact", "fiber_diam", "type", "inner", "sample"])
analysis = grouped.agg({"threshold": [np.mean, rangecalc]})
analysis.columns = ["_".join(col_name).rstrip("_") for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)
# normalize all range values to the max for that fiber diameter
analysis['threshold_range_norm'] = analysis.groupby('fiber_diam')['threshold_rangecalc'].transform(
    lambda x: x / x.max()
)


plt.figure()
g = sns.pointplot(
    data=analysis,
    y="threshold_range_norm",
    x="fiber_diam",
    hue="type",
    palette=pal2d3d,
    dodge=True,
)
plt.title("intrafascicle")
plt.ylabel("Normalized Threshold Range")
plt.xlabel("Fiber Diameter (μm)")
# plt.yscale('log')
plt.gca().get_legend().set_title("")
plt.gcf().set_size_inches([3, 2])
plt.ylim(0, None)

# significance using wilcoxon
meltmatch = datamatch_merge(
    analysis.query('type=="extrusion"'),
    analysis.query('type=="true-3D"'),
    "threshold_range_norm",
    merge_cols=["contact", "fiber_diam", "inner", "sample"],
).drop(columns="type")
from scipy.stats import wilcoxon

for diam in meltmatch.fiber_diam.unique():
    extrusion = meltmatch.query(f'fiber_diam=={diam}')['threshold_range_norm']
    true3d = meltmatch.query(f'fiber_diam=={diam}')['threshold_range_norm3d']
    print(f'{diam} μm: {wilcoxon(extrusion,true3d)}')


# now do variance between fascicle mean thresholds
grouped = analysis.groupby(["contact", "fiber_diam", "type", "sample"])
analysis = grouped.agg({"threshold_mean": [rangecalc]})
analysis.columns = ["_".join(col_name).rstrip("_") for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)
# normalize all range values to the max for that fiber diameter
analysis['threshold_mean_range_norm'] = analysis.groupby('fiber_diam')['threshold_mean_rangecalc'].transform(
    lambda x: x / x.max()
)
plt.figure()
g = sns.pointplot(
    data=analysis,
    y="threshold_mean_range_norm",
    x="fiber_diam",
    hue="type",
    palette=pal2d3d,
    dodge=True,
    estimator=estimator,
    errorbar=errorbar,
)
plt.title("interfascicle")
plt.ylabel("Normalized Threshold Range")
plt.xlabel("Fiber Diameter (μm)")
plt.gca().get_legend().set_title("")
plt.ylim(0, None)
plt.gcf().set_size_inches([3, 2])

# threshold variances for whole nerve
vardat = repeated.copy()
grouped = vardat.groupby(["contact", "fiber_diam", "sim", "type", "sample"])
analysis = grouped.agg({"threshold": [rangecalc]})
analysis.columns = ["_".join(col_name).rstrip("_") for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)
# normalize all range values to the max for that fiber diameter
analysis['threshold_range_norm'] = analysis.groupby('fiber_diam')['threshold_rangecalc'].transform(
    lambda x: x / x.max()
)

plt.figure()
g = sns.pointplot(
    data=analysis,
    y="threshold_range_norm",
    x="fiber_diam",
    hue="type",
    palette=pal2d3d,
    dodge=True,
    # legend=True
)
plt.title("intra-sample")
plt.ylabel("Normalized Threshold Range")
plt.xlabel("Fiber Diameter (μm)")
plt.gca().get_legend().set_title("")
plt.ylim(0, None)
plt.gcf().set_size_inches([3, 2])
# %% CCC comparison maincomp
sns.set(context="paper", style="whitegrid")
diam = 3
comp = cath_comparison
deformation = "Structural"
nsimdata = concats.query(
    f"fiber_diam in [3] and contact in {comp} and deformation==@deformation"
)  # TODO replace all cath comparison with non
g = sns.relplot(
    data=nsimdata.rename(columns={"nerve_label": "Sample"}),
    kind="scatter",
    x="threshold",
    y="threshold3d",
    # hue='Sample',
    color="white",
    s=10,
    facet_kws={"sharex": False, "sharey": False, "margin_titles": True},
    edgecolor="black",
    linewidth=0.5,
    alpha=1,
    label="single fiber",
)
lim = {}
ax = g.axes[0, 0]
# TODO Clean up this calc
rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
print(r)
perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
lim[deformation] = np.amax([ax.get_xlim()[1], ax.get_ylim()[1]])
# add correlation to plot
# ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
print(f"{diam} {deformation} {comp[0]} μm CCC: {r ** 2:.2f}")
# ax.set_title(f'{diam} μm')
ax.plot(
    [0, lim[deformation]],
    [0, lim[deformation]],
    "--k",
    linewidth=1,
    # label="unity line",
)
plt.legend(bbox_to_anchor=(0.9, 1.25))
# ax.set_aspect('equal', 'box')
# ax.apply_aspect()
ax.set_xlim([0, lim[deformation]])
ax.set_ylim([0, lim[deformation]])

# ax.set_yticks(ax.get_xticks())
g.set_titles("D: {col_name} μm")
g.set_xlabels("extrusion threshold (mA)")
g.set_ylabels("true-3D threshold (mA)")
g.set_titles(row_template="", col_template="Deformation: {col_name}")

g.set_titles(row_template="", col_template="Deformation: {col_name}")
mid = [np.diff(plt.xlim()) / 2, np.diff(plt.ylim()) / 2]
mid = [float(x) for x in mid]
plt.arrow(mid[0] - 0.25, mid[1] + 0.25, -0.75, 0.75, color="black", width=0.08)
plt.text(mid[0] - 1.8, mid[1] + 1.4, "true-3D higher")
plt.arrow(mid[0] + 0.25, mid[1] - 0.25, +0.75, -0.75, color="black", width=0.08)
plt.text(mid[0], mid[1] - 1.7, "extrusion higher")
plt.gcf().set_size_inches(1.5, 1.5)
# ax.apply_aspect()
ax.set_xlim([0, lim[deformation]])
ax.set_ylim([0, lim[deformation]])
plt.xticks([0, 1, 2, 3])
# %% dose-response compare
peses = []
pemeans = []
onsets_sats = {}
for comparison in comparisons:
    for stringdat in ["Undeformed", "Structural", "ASCENT"]:
        thiscontact = comparison[0]
        subdat = newdefdat.query(f"deformation=='{stringdat}' and contact in {comparison}")
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
        assert thiscontact != np.nan
        peses.append(pes)
        print("Max 3 um", stringdat, np.amax(pes.query("fiber_diam==3").pe))
        print("Max 13 um", stringdat, np.amax(pes.query("fiber_diam==13").pe))

        # now calculate percent error for population onset and saturation
        pemean = []
        for level in ["onset", "half", "saturation"]:
            for fiber_diam in compiled_data.fiber_diam.unique():
                a = compiled_data.query(f"level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam}")[
                    "threshold"
                ].values
                b = compiled_data.query(f"level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}")[
                    "threshold"
                ].values
                pe_res = pe_noabs(np.median(b), np.median(a), doabs=False)
                pemean.append({"level": level, "fiber_diam": fiber_diam, "pe": pe_res})

        pemean = pd.DataFrame(pemean)
        pemean["deformation"] = stringdat
        pemean["contact"] = thiscontact
        pemeans.append(pemean)
        print("Max 3 um median", stringdat, np.amax(pemean.query("fiber_diam==3").pe))
        print("Max 13 um median", stringdat, np.amax(pemean.query("fiber_diam==13").pe))

        onsets_sats[stringdat] = compiled_data

sns.set(style="whitegrid", context="paper")
allpes = pd.concat(peses).query('deformation != "Undeformed"')
# allpes['level'].replace({
#     'onset':'10','half':'50','saturation':'90'},inplace=True)
# allpes['deformation'].replace({
#     'ASCENT':'ASCENT','Structural':'Structural','Undeformed':'Undef.'},inplace=True)
# allpes['comb'] = allpes['level'].astype(str) + '\n' + allpes['deformation']
# allpes.sort_values(by=['deformation','level'])
allpes["level"] = pd.Categorical(allpes.level, categories=["onset", "half", "saturation"], ordered=True)
allpes["deformation"] = pd.Categorical(allpes.deformation, categories=["Structural", "ASCENT"], ordered=True)
allpes["contact"] = pd.Categorical(allpes.contact, categories=["anodic", "center", "cathodic"], ordered=True)
allpes.sort_values(by=["contact", "deformation"], inplace=True)
allpes["comb"] = allpes["deformation"].astype(str) + "\n" + allpes["contact"].astype(str)

sns.set(context='paper', style='white', font_scale=1)
# allpes=allpes.query('deformation!="ASCENT"')
# allpes["deformation"] = allpes.deformation.replace({'Structural':'Deformed'})
# allpes["deformation"] = pd.Categorical(
#     allpes.deformation, categories=["Undeformed", "Deformed"], ordered=True
# )
sns.set(context='paper', style='white', font_scale=1)
sns.catplot(
    data=allpes.query(f'contact in {cath_comparison}'),
    kind="strip",
    row="deformation",
    y="pe",
    x="level",
    hue="fiber_diam",
    # col="nerve_label",
    palette='RdPu',
    margin_titles=True,
    dodge=True,
    linewidth=1,
)
for ax in plt.gcf().axes:
    ax.set_xlabel("")
    ax.axhline(0, ls="--", color="black")

# redo with facetgrid that has barplot and stripplot
plt.figure()
g = sns.FacetGrid(
    data=allpes.query(f'contact in {cath_comparison}'),
    row="deformation",
    # col="nerve_label",
    palette='RdPu',
    margin_titles=True,
)
for ax in plt.gcf().axes:
    ax.set_xlabel("")
    ax.axhline(0, ls="-", color="black")
g.map_dataframe(sns.stripplot, y="pe", x="level", dodge=True, linewidth=1, hue='fiber_diam', palette='RdPu')
g.map_dataframe(sns.barplot, y="pe", x="level", hue='fiber_diam', facecolor='white', edgecolor='k', errorbar=None)
g.set_ylabels("Percent Error")
plt.gcf().set_size_inches(4, 8 / 3)
rownames(g, row_template="{row_name}")
g.axes[1, 0].set_ylabel(f'Relative Percent Difference (%)\n{g.axes[1,0].get_ylabel()}')
g.set_titles(row_template='')
for ax in g.axes.ravel():
    ax.tick_params(pad=0)
    ax.set_yticks([-40, -20, 0, 20, 40])
plt.ylim(-50, 50)
g.set_xticklabels(['onset (10%)', 'half (50%)', 'saturation (90%)'])
# %% now for only cathodic calculate the median and plot barplot with values
specdat = allpes.query(f'contact in {cath_comparison}')
specdat = specdat.groupby(['deformation', 'level', 'fiber_diam']).agg({'pe': 'median'}).reset_index()
g = sns.catplot(
    row='deformation', col='level', x='fiber_diam', y='pe', data=specdat, kind='bar', palette='RdPu', margin_titles=True
)
# with values
for ax in g.axes.ravel():
    for i in ax.containers:
        ax.bar_label(i, fmt='%.1f')
plt.gcf().set_size_inches(4, 4)
# %%
sns.set(style='white')
for samplenum, prefix in zip([2521, 3721, 5721, 6721], ['2L', '3R', '5R', '6R']):
    # for samplenum,prefix in zip([3721,],['3R']):
    criteria = {
        'partial_matches': True,
        'include_downstream': False,
        'indices': {'sample': [samplenum], 'model': None, 'sim': None},
    }

    q = Query(criteria)
    q.run()

    results = q.summary()

    sample_index = results['samples'][0]['index']

    fig, ax = plt.subplots(1, 1)
    item: Sample = q.get_object(Object.SAMPLE, [results['samples'][0]['index']])
    slide = item.slides[0]
    cmap = cm.get_cmap('cmr.gem', len(slide.inners())).reversed()
    # plot fascicle ECD

    ECD = [{'inner': i, "D": f.ecd()} for i, f in enumerate(slide.inners())]
    ECD = pd.DataFrame(ECD)
    ECD.sort_values('D', inplace=True)
    colors = [to_hex(cmap(i)) for i in range(len(slide.inners()))]
    ECD['color'] = None
    for i, row in enumerate(ECD.itertuples()):
        ECD.loc[row.Index, 'color'] = colors[i]
    ECD.sort_values('inner', inplace=True)

    slide.plot(
        fix_aspect_ratio=True,
        final=False,
        ax=ax,
        inner_index_labels=True,
        # scalebar=True,
        # scalebar_length=100,
        # scalebar_units='μm',
        fascicle_colors=ECD['color'].values,
        inner_format="w-",
        line_kws=dict(linewidth=1),
        colors_for_outers=True,
        inners_flag=True,
    )
    plt.xlabel('\u03bcm')
    plt.ylabel('\u03bcm')

    # load contact coordinates
    dire = rf'D:\threed_final\input\contact_coords\{prefix}defDS5'
    for i in [1]:  # NOTE set this to only load caudal coords
        contact_coords = np.loadtxt(dire + r'\pcs' + str(i + 1) + '.txt', skiprows=8)[:, :2]
        plt.scatter(contact_coords[:, 0], contact_coords[:, 1], color='k', label='contacts' if i == 0 else '_')
        plt.legend()
    plt.show()
    fname = str(sample_index)
    fmt = 'svg'

    # plot sxcatterplot of ECD with x=0 for all
    fig, ax = plt.subplots(1, 1)
    for row in ECD.itertuples():
        ax.scatter(0, row.D, color=row.color, marker='s', s=100)
        ax.scatter(0, row.D, color='black', marker='_', s=100)
    plt.gca().set_aspect(0.0015)
    plt.xticks([], [])
    plt.ylabel('fascicle diameter (μm)')
    plt.ylim(0, 1000)
    plt.show()

    sns.set(context='paper', style='white', font_scale=1)
    plt.figure()
    drthis = defdr.query(
        f"fiber_diam in [3] and contact in {cath_comparison} and nerve_label ==@prefix and deformation=='Structural'"
    )

    # Query the rows with type '3D'
    df_3d = drthis.query("modeltype == 'true-3D'")

    # Query the rows with type '2D' and sample ending in '1'
    df_2d = drthis.query("modeltype == 'extrusion'")

    # Merge the 3D and 2D data, keeping track of original row indices
    merged_df = pd.merge(
        df_3d,
        df_2d,
        on=["nerve_label", "master_fiber_index", "nsim", "deformation"],
        suffixes=("_3d", "_2d"),
        how="left",
    )  # TODO remove this how

    # Update the 'inner', 'outer', and 'fiber' columns in the original DataFrame
    drthis.loc[df_3d.index, "inner"] = merged_df["inner_2d"].values
    drthis.loc[df_3d.index, "outer"] = merged_df["outer_2d"].values
    drthis.loc[df_3d.index, "fiber"] = merged_df["fiber_2d"].values

    g = sns.scatterplot(
        data=drthis,
        y="percent_activated",
        x="threshold",
        hue="inner",
        palette=ECD['color'].values,
        linewidth=0,
        s=50,
        marker=r"$\backslash$",
        # legend=False
        hue_order=range(len(slide.fascicles)),
    )

    sns.lineplot(
        data=drthis,
        y="percent_activated",
        x="threshold",
        color="k",
        linewidth=2,
        estimator=None,
        style="modeltype",
        legend=False,
        zorder=1,
        alpha=1,
    )

    plt.xlabel("Threshold (mA)")
    plt.ylim([0, 1])
    # plt.gcf().set_size_inches(8,4)
    plt.ylabel("Proportion Activated")
    # create legend, circle = extrusion, X= true-3D
    # create handles manually
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], label="extrusion", color="k", linestyle="-", alpha=0.6),
        Line2D([0], [0], label="true-3D", color="k", linestyle="--", alpha=0.6),
    ]
    legend_labels = ["extrusion", "true-3D"]
    g.legend(handles=legend_elements, labels=legend_labels, loc="lower right")
    plt.gcf().set_size_inches(2, 2)

    plt.figure()
    g = sns.swarmplot(
        data=drthis,
        y="threshold",
        x="modeltype",
        hue="inner",
        palette=ECD['color'].values,
        linewidth=0,
        s=2,
        hue_order=range(len(slide.fascicles)),
        # legend=False
    )
    plt.xlabel("")
    plt.ylabel("Threshold (mA)")
    plt.legend([], [], frameon=False)
    plt.ylim(0, None)
    thismatch = datamatch_merge(
        drthis.query('modeltype=="extrusion"'),
        drthis.query('modeltype=="true-3D"'),
        "threshold",
        merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
    ).drop(columns="modeltype")
    plt.gcf().set_size_inches(2, 2)

    # add fascicle size to drthis. first turn ECD into a dict mapping between inner and D
    ECDdict = ECD.set_index('inner')['D'].to_dict()
    # now add to drthis
    drthis['inner_ecd'] = drthis['inner'].replace(ECDdict)

    # now plot
    plt.figure()
    sns.lineplot(
        data=drthis,
        x="inner_ecd",
        y="threshold",
        hue="modeltype",
        palette=pal2d3d,
        # linewidth=1,
        # marker='o'
    )

    plt.figure()
    sns.scatterplot(
        data=thismatch,
        x="threshold",
        y="threshold3d",
        hue="inner",
        palette=ECD['color'].values,
        legend=False,
        hue_order=range(len(slide.fascicles)),
        s=10,
    )
    plt.xlim(0, 4)
    plt.ylim(0, 4)
    plt.gca().set_aspect("equal")
    plt.plot([0, 5], [0, 5], "k--")
    plt.ylabel("True-3D Threshold")
    plt.xlabel("Extrusion threshold")
    plt.gcf().set_size_inches(2, 2)
# %% threshold variances and coefficient of variation intrafascicle and inter
estimator, errorbar = "median", ('ci', 95)
from scipy.stats import variation

sns.set(font_scale=1, style="white", context='paper')
vardat = repeated_deformation.query('contact=="cathodic" and deformation=="Structural"')
grouped = vardat.groupby(["contact", "fiber_diam", "type", "inner", "sample"])
analysis = grouped.agg({"threshold": [np.var, np.mean, variation]})
analysis.columns = ["_".join(col_name).rstrip("_") for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)


plt.figure()
g = sns.pointplot(
    data=analysis,
    y="threshold_variation",
    x="fiber_diam",
    hue="type",
    palette=pal2d3d,
    dodge=True,
    estimator=estimator,
    errorbar=errorbar,
)
plt.title("intrafascicle")
plt.ylabel("Threshold CoV")
plt.xlabel("Fiber Diameter (μm)")
# plt.yscale('log')
plt.gca().get_legend().set_title("")
plt.gcf().set_size_inches([6, 4])
plt.ylim(0, None)
plt.gcf().set_size_inches(3, 2)

# test significance non-parametric
from scipy.stats import mannwhitneyu

for diam in analysis.fiber_diam.unique():
    extrusion = analysis.query(f'fiber_diam == {diam} and type == "extrusion"').threshold_variation
    true3d = analysis.query(f'fiber_diam == {diam} and type == "true-3D"').threshold_variation
    print(f'{diam} μm')
    print(mannwhitneyu(extrusion, true3d))

# now do variance between fascicle mean thresholds
grouped = analysis.groupby(["contact", "fiber_diam", "type", "sample"])
analysis = grouped.agg({"threshold_mean": [np.var, variation, np.count_nonzero]})
analysis.columns = ["_".join(col_name).rstrip("_") for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)

plt.figure()
g = sns.pointplot(
    data=analysis,
    y="threshold_mean_variation",
    x="fiber_diam",
    hue="type",
    palette=pal2d3d,
    dodge=True,
    estimator=estimator,
    errorbar=errorbar,
)
plt.title("interfascicle")
plt.ylabel("Threshold CoV")
plt.xlabel("Fiber Diameter (μm)")
plt.gca().get_legend().set_title("")
plt.ylim(0, None)
plt.gcf().set_size_inches(3, 2)

# significance
for diam in analysis.fiber_diam.unique():
    extrusion = analysis.query(f'fiber_diam == {diam} and type == "extrusion"').threshold_mean_variation
    true3d = analysis.query(f'fiber_diam == {diam} and type == "true-3D"').threshold_mean_variation
    print(f'{diam} μm')
    print(mannwhitneyu(extrusion, true3d))

# threshold variances for whole nerve
vardat = repeated_deformation.query('contact=="cathodic" and deformation=="Structural"')
grouped = vardat.groupby(["contact", "fiber_diam", "sim", "type", "sample"])
analysis = grouped.agg({"threshold": [np.var, np.mean, variation]})
analysis.columns = ["_".join(col_name).rstrip("_") for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)

plt.figure()
g = sns.pointplot(
    data=analysis,
    y="threshold_variation",
    x="fiber_diam",
    hue="type",
    palette=pal2d3d,
    dodge=True,
    estimator=estimator,
    errorbar=errorbar,
)
plt.title("intra-sample")
plt.ylabel("Threshold CoV")
plt.xlabel("Fiber Diameter (μm)")
plt.gca().get_legend().set_title("")
plt.ylim(0, None)
plt.gcf().set_size_inches(3, 2)
# %% dose-response
sns.set(context='paper', style='ticks')
for nerve in ["2L", '3R', '5R', '6R']:
    alldr = defdr.query(f"nerve_label in {defsamples}").query(f"nerve_label=='{nerve}'")
    alldr["deformed"] = alldr["deformation"] != "Undeformed"
    allstore = []
    for contact in alldr.contact.unique():
        if contact == "3D":
            continue
        thisdat = alldr.query('contact in [@contact, "3D"]')
        thisdat["contact"] = contact
        allstore.append(thisdat)
    alldr = pd.concat(allstore)
    alldr["deformation"] = pd.Categorical(
        alldr["deformation"],
        ordered=True,
        categories=["Undeformed", "Structural", "ASCENT"],
    )

    plt.figure()
    g = sns.relplot(
        kind="line",
        style="deformation",
        data=alldr.query("fiber_diam in [3]"),
        y="percent_activated",
        x="threshold",
        units="nerve_label",
        hue="modeltype",
        palette=pal2d3d,
        estimator=None,
        linewidth=2,
        facet_kws={"sharex": True, "margin_titles": True},
        row="deformed",
        col="contact",
        legend=True,
    )

    g.legend.set_title("")

    for ax in g.axes.ravel():
        ax.set_ylim([0, 1])
        ax.set_xlim([0, None])
    g.set_ylabels("Proportion fibers active")
    g.set_titles(col_template="{col_name}", row_template='')
    plt.subplots_adjust(hspace=0.15)

    g.set_xlabels("Threshold (mA)")
    # change the line width for the legend
    for line, l in zip(g.legend.get_lines(), g.legend.get_texts()):
        line.set_linewidth(2.0)
        if l.get_text() in ['deformation', 'modeltype']:
            l.set_text('')
    for ax in g.axes.ravel():
        for loc in [0.1, 0.5, 0.9]:
            ax.axhline(loc, color="gray", linestyle="-", alpha=0.5, linewidth=1)
    sns.move_legend(g, [0.75, 0.2], facecolor='white', framealpha=1, frameon=True, edgecolor='black')
    g.axes[0][0].set_ylabel("Active fibers (%)\nUndeformed")
    g.axes[1][0].set_ylabel("Active fibers (%)\nDeformed")
    # rownames(g, row_template="{row_name}-slice\nProportion fibers active")
    g.axes[0][0].set_yticks([0, 0.1, 0.5, 0.9, 1], ['0', '10', '50', '90', '100'])
    g.fig.set_size_inches(4, 2.5)
    plt.suptitle(nerve, y=1.1)
# %% plot activation order
sns.set(font_scale=1, context='paper', style="white")
for nerve in ["2L", '3R', '5R', '6R']:
    plotactidata = newdefdr.query(f"fiber_diam in [3] and nerve_label=='{nerve}'")
    plt.figure()
    g = sns.catplot(
        kind="swarm",
        row="fiber_diam",
        data=plotactidata,
        y="percent_activated",
        x="contact",
        units="nerve_label",
        palette="plasma",
        hue="percent_activated3d",
        # hue='inner',
        estimator=None,
        linewidth=0,
        # facet_kws={"margin_titles": True},
        s=10,
    )
    plt.subplots_adjust(top=0.87)
    # plt.suptitle(stringdat, x=0.37)
    g.set_titles(row_template="", col_template="{col_name}")
    g.axes[0][0].set_xlabel("")
    g.axes[0][0].set_ylabel("Proportion fibers activated")
    g.set_xlabels("")
    g.legend.remove()
    norm = plt.Normalize(0, 1)
    sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
    sm.set_array([])
    plt.axvline(0.5, linestyle="--", color="black", alpha=1)
    # Remove the legend and add a colorbar
    g.figure.colorbar(
        sm,
        ax=g.axes.ravel().tolist(),
        aspect=10,
        shrink=0.8,
        label="Proportion true-3D\nfibers activated",
        pad=0.1,
    ).ax.yaxis.set_ticks_position("left")
    g.fig.set_size_inches(14, 4)
    for i, con in enumerate(newdefdr.contact.sort_values().unique()[1:]):
        shortdat = plotactidata.query("contact==@con")
        data2d = shortdat.sort_values("percent_activated").master_fiber_index
        data3d = shortdat.sort_values("percent_activated3d").master_fiber_index
        rc = compute_reorder_cost(list(data2d), list(data3d))
        g.axes[0][0].text(i + 0.6, 1.065, f'AR: {round(rc, 3)}')
    plt.gcf().set_size_inches(6, 2)
    # now, calculate which fibers have a change in activation order between 2D and 3D, plot a grey line between the new and old values
    for ax in g.axes.ravel():
        for i, con in enumerate(newdefdr.contact.sort_values().unique()[1:]):
            shortdat = plotactidata.query("contact==@con")
            for fiber_data in shortdat.itertuples():
                if fiber_data.percent_activated != fiber_data.percent_activated3d:
                    ax.plot(
                        [i + 0.55, i + 1],
                        [fiber_data.percent_activated, fiber_data.percent_activated3d],
                        color="gray",
                        alpha=0.5,
                        linewidth=0.5,
                    )
    # sys.exit()
# %% activation order fiberdiam
newdefdr = defdr.copy()
newdefdr = newdefdr.query("'6R' in nerve_label and deformation!='ASCENT'").sort_values("modeltype")
newdefdr["deformation"] = pd.Categorical(newdefdr["deformation"], categories=["Undeformed", "Structural"])
sns.set(font_scale=1.5)
sns.set_style("whitegrid")
newdefdr["percent_activated3"] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newdefdr.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        f"fiber_diam == 3 and nerve_label == @row.nerve_label and modeltype == @row.modeltype and sim == @row.sim and master_fiber_index == @row.master_fiber_index and contact in {cath_comparison} and deformation==@row.deformation"
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert val is not np.nan
    newdefdr.loc[row.Index, "percent_activated3"] = val
#  plot
sns.set(font_scale=1, style="white", context='paper')
plt.figure()
newdefdr['fiber_diam'] = newdefdr['fiber_diam'].astype(str)
g = sns.catplot(
    kind="swarm",
    # row='deformation',
    data=newdefdr.query(f"fiber_diam in ['3','13'] and contact in {cath_comparison} and deformation=='Structural'"),
    y="percent_activated",
    x="fiber_diam",
    units="nerve_label",
    col="modeltype",
    palette="plasma",
    hue="percent_activated3",
    # hue='inner',
    estimator=None,
    linewidth=0,
    margin_titles=True,
    s=10,
)
xlim = plt.xlim()
# g.map_dataframe(sns.lineplot,
#                 y="percent_activated",
#                 x="fiber_diam",
#                 # hue='inner',
#                 estimator=None,
#                 units = 'master_fiber_index',
#                 color='black',
#                 alpha=0.3
#                 )
newdefdr['fiber_diam'] = newdefdr['fiber_diam'].astype(int)

plt.subplots_adjust(top=0.87)
# plt.suptitle(stringdat, x=0.37)
g.set_titles(row_template="", col_template="{col_name}")
g.axes[0][0].set_xlabel("")
g.axes[0][0].set_ylabel("Proportion fibers activated")
g.set_xlabels("D (μm)")
g.legend.remove()
norm = plt.Normalize(0, 1)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])

# Remove the legend and add a colorbar
g.figure.colorbar(
    sm,
    ax=g.axes.ravel().tolist(),
    aspect=10,
    shrink=0.8,
    label="Proportion 3 μm\nfibers activated",
    pad=0.1,
).ax.yaxis.set_ticks_position("left")

for i, mod in enumerate(newdefdr.modeltype.sort_values().unique()):
    shortdat = newdefdr.query(
        f"modeltype==@mod and fiber_diam in [13] and contact in {cath_comparison} and deformation=='Structural'"
    )
    data2d = shortdat.sort_values("percent_activated").master_fiber_index
    data3d = shortdat.sort_values("percent_activated3").master_fiber_index
    rc = compute_reorder_cost(list(data2d), list(data3d))
    ax = g.axes[0][i]
    text = ax.get_title() + f" - AR: {round(rc,3)}"
    ax.set_title(text)
plt.gcf().set_size_inches([5, 2])
plt.xlim(xlim)

# now do the same for this one (gray lines between changed points)
for ax in g.axes.ravel():
    for i, mod in enumerate(newdefdr.modeltype.sort_values().unique()):
        shortdat = newdefdr.query(
            f"modeltype==@mod and fiber_diam in [13] and contact in {cath_comparison} and deformation=='Structural'"
        )
        for fiber_data in shortdat.itertuples():
            if fiber_data.percent_activated != fiber_data.percent_activated3:
                ax.plot(
                    [0.55, 1],
                    [fiber_data.percent_activated, fiber_data.percent_activated3],
                    color="gray",
                    alpha=0.5,
                    linewidth=0.5,
                )
# %% dose-response  each
sns.set(context='paper', style='ticks')
alldr = defdr.query(f"nerve_label in {defsamples} and deformation=='Structural' and contact in {cath_comparison}")

plt.figure()
g = sns.relplot(
    kind="line",
    # style="deformation",
    data=alldr.query("fiber_diam in [3,13]"),
    y="percent_activated",
    x="threshold",
    units="nerve_label",
    hue="modeltype",
    palette=pal2d3d,
    estimator=None,
    # linewidth=2,
    facet_kws={"sharex": False, "margin_titles": True},
    row="deformed",
    col="fiber_diam",
    alpha=0.5,
)

g.legend.set_title("")

for ax in g.axes.ravel():
    ax.set_ylim([0, 1])
    ax.set_xlim([0, None])
g.set_ylabels("Proportion fibers active")
g.set_titles(col_template="{col_name}", row_template='')
plt.subplots_adjust(hspace=0.15)

g.set_xlabels("Threshold (mA)")
# change the line width for the legend
for line, l in zip(g.legend.get_lines(), g.legend.get_texts()):
    line.set_linewidth(2.0)
    # if l.get_text() in ['deformation', 'modeltype']:
    #     l.set_text('')
for ax in g.axes.ravel():
    for loc in [0.1, 0.5, 0.9]:
        ax.axhline(loc, color="gray", linestyle="-", alpha=0.5, linewidth=1)
sns.move_legend(g, [0.45, 0.3], facecolor='white', framealpha=1, frameon=True, edgecolor='black')
g.axes[0][0].set_ylabel("Active fibers (%)")
g.set_titles(col_template='D: {col_name} μm', row_template='')
# g.axes[1][0].set_ylabel("Active fibers (%)\nDeformed")
# rownames(g, row_template="{row_name}-slice\nProportion fibers active")
g.axes[0][0].set_yticks([0, 0.1, 0.5, 0.9, 1], ['0', '10', '50', '90', '100'])

g.fig.set_size_inches(4, 2)

# calculate median for each diameter and modeltype, plot on each subplot
meddata = []
data_bynerve = []
for activation in np.linspace(0, 1, 50):
    # for each activation level, calculate the median threshold for each modeltype and fiber diameter
    for modeltype in alldr.modeltype.unique():
        for fiber_diam in alldr.fiber_diam.unique():
            # take the highest threshold value that is less than or equal to the activation level, then take the median of these
            thisdat = alldr.query(
                "modeltype==@modeltype and fiber_diam==@fiber_diam and percent_activated<=@activation"
            )
            # take the max threshold for each nerve
            thisdat = thisdat.groupby(['nerve_label']).agg({'threshold': 'max'}).reset_index()
            data_bynerve.append(thisdat)
            data_bynerve[-1]['fiber_diam'] = fiber_diam
            data_bynerve[-1]['modeltype'] = modeltype
            meddata.append(
                {
                    'modeltype': modeltype,
                    'fiber_diam': fiber_diam,
                    'percent_activated': activation,
                    'threshold': thisdat.median().threshold,
                }
            )
meddata = pd.DataFrame(meddata).dropna()
for ax, diam in zip(g.axes.ravel(), [3, 13]):
    for modeltype, color in zip(['extrusion', 'true-3D'], pal2d3d):
        thisdat = meddata.query(f"fiber_diam=={diam} and modeltype==@modeltype")
        sns.lineplot(
            data=thisdat,
            y='percent_activated',
            x='threshold',
            ax=ax,
            color=color,
            estimator=None,
            linewidth=2,
            legend=False,
        )
data_bynerve = pd.concat(data_bynerve).dropna()

# plot the error from true-3D median to extrusion median and for each sample as well
matched = datamatch_merge(
    meddata.query('modeltype=="extrusion"'),
    meddata.query('modeltype=="true-3D"'),
    'threshold',
    merge_cols=['fiber_diam', 'percent_activated'],
).dropna()
# calculate the percent error
matched['pe'] = pe_noabs(matched.threshold3d, matched.threshold, doabs=False)
plt.figure()
g = sns.relplot(
    data=matched.query('fiber_diam in [3,13]'),
    kind='line',
    x='percent_activated',
    y='pe',
    hue='fiber_diam',
    color='black',
)
# add to plot
for ax, diam in zip(g.axes.ravel(), [3, 13]):
    for modeltype, color in zip(['extrusion', 'true-3D'], pal2d3d):
        thisdat = data_bynerve.query(f"fiber_diam=={diam} and modeltype==@modeltype")
        sns.lineplot(
            data=thisdat, y='threshold', x='nerve_label', ax=ax, color=color, estimator=None, linewidth=1, legend=False
        )
# %% plot error
meddata = []
# redo with just onset half and sat, and plot with barplot
for activation in [0.1, 0.5, 0.9]:
    # for each activation level, calculate the median threshold for each modeltype and fiber diameter
    for modeltype in alldr.modeltype.unique():
        for fiber_diam in alldr.fiber_diam.unique():
            # take the highest threshold value that is less than or equal to the activation level, then take the median of these
            thisdat = alldr.query(
                "modeltype==@modeltype and fiber_diam==@fiber_diam and percent_activated<=@activation"
            )
            # take the max threshold for each nerve
            thisdat = thisdat.groupby(['nerve_label']).agg({'threshold': 'max'}).reset_index()
            meddata.append(
                {
                    'modeltype': modeltype,
                    'fiber_diam': fiber_diam,
                    'percent_activated': activation,
                    'threshold': thisdat.median().threshold,
                }
            )
# plot
meddata = pd.DataFrame(meddata).dropna()
plt.figure()
matched = datamatch_merge(
    meddata.query('modeltype=="extrusion"'),
    meddata.query('modeltype=="true-3D"'),
    'threshold',
    merge_cols=['fiber_diam', 'percent_activated'],
).dropna()
# calculate the percent error
matched['pe'] = pe_noabs(matched.threshold3d, matched.threshold, doabs=False)
sns.barplot(
    data=matched,
    x='percent_activated',
    y='pe',
    hue='fiber_diam',
    palette='RdPu'
    # black outline
    ,
    edgecolor='black',
    linewidth=1,
)
plt.legend(ncol=2, title='D (μm)')
plt.axhline(0, ls='-', color='black')
plt.xlabel('')
plt.xticks([0, 1, 2], ['onset (10%)', 'half (50%)', 'saturation (90%)'])
plt.ylabel('Relative Percent Difference (%)')
plt.gcf().set_size_inches(4, 3)

# %% mct fiberdiam order
sns.set(font_scale=1, context='paper', style="whitegrid")
scores = []
for cc in sorted(imdata["active_src_index"].unique()):
    shortdat = imdata.query(f'active_src_index=="{cc}"')
    for nerve_label in pd.unique(shortdat["nerve_label"]):
        for modeltype in pd.unique(shortdat["type"]):
            finallythedat = shortdat.query(f'nerve_label=="{nerve_label}" and type=="{modeltype}"')
            datasmol = finallythedat.query('fiber_diam==3').sort_values('threshold').master_fiber_index
            databeeg = finallythedat.query('fiber_diam==13').sort_values('threshold').master_fiber_index
            rc = compute_reorder_cost(list(datasmol), list(databeeg))
            scores.append(
                {
                    "active_src_index": cc,
                    "score2d3d": rc,
                    "nerve_label": nerve_label,
                    "modeltype": modeltype,
                }
            )
scoredat = pd.DataFrame(scores)
scoredatnew = scoredat.copy()
scoredatnew["active_src_index"] = scoredatnew["active_src_index"].replace(
    {"1": "MCT", "2": "MCT", "3": "MCT", "0": "MCT"}
)
scoredatnew["active_src_index"] = scoredatnew["active_src_index"].replace(
    {"MCT": "Multi-Contact Cuff", "LN": "LivaNova Cuff"}
)
plt.figure()
g = sns.barplot(
    data=scoredatnew,
    y="score2d3d",
    x="active_src_index",
    hue="modeltype",
    palette=pal2d3d,
    # edgecolor='k',
    errorbar=('sd', 1),
)
sns.stripplot(
    data=scoredatnew,
    y="score2d3d",
    x="active_src_index",
    hue='modeltype',
    dodge=True,
    edgecolor='black',
    linewidth=0.5,
    s=3,
    palette=pal2d3d,
    jitter=0.25,
)

plt.xlabel("")
plt.ylabel("Activation Reordering")
plt.legend(title='', ncol=2)
plt.gcf().set_size_inches(2.5, 1.5)
plt.ylim(0, 0.6)
plt.axhline(0.45, color='red', ls='--')
plt.text(-0.3, 0.5, 'random', color='red')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title="", handles=handles[:2], labels=labs[:2], ncol=1, loc='upper right')
# %%# %% plot activation order newmethod
sns.set(font_scale=1, context='paper', style="white")
for nerve in ["2L", '3R', '5R', '6R']:
    plotactidata = newdefdr.query(f"fiber_diam in [3] and nerve_label=='{nerve}'")
    plt.figure()
    g = sns.relplot(
        kind="scatter",
        row="fiber_diam",
        data=plotactidata,
        x="percent_activated",
        col="contact",
        units="nerve_label",
        palette="plasma",
        y="percent_activated3d",
        hue='percent_activated3d',
        # hue='inner',
        # estimator=None,
        # linewidth=0,
        # facet_kws={"margin_titles": True},
        s=10,
    )
    plt.subplots_adjust(top=0.87)
    # plt.suptitle(stringdat, x=0.37)
    g.set_titles(row_template="", col_template="{col_name}")
    g.axes[0][0].set_xlabel("")
    g.axes[0][0].set_ylabel("Proportion fibers activated")
    g.set_xlabels("")
    # g.legend.remove()
    norm = plt.Normalize(0, 1)
    sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
    sm.set_array([])

    for i, con in enumerate(newdefdr.contact.sort_values().unique()[1:]):
        shortdat = plotactidata.query("contact==@con")
        data2d = shortdat.sort_values("percent_activated").master_fiber_index
        data3d = shortdat.sort_values("percent_activated3d").master_fiber_index
        rc = compute_reorder_cost(list(data2d), list(data3d))
        g.axes[0][i + 1].set_title(f'AR: {round(rc, 3)}\ncomparison: {con}')
    plt.gcf().set_size_inches(8, 2)
    g.set_xlabels('Percent activated comparison')
    g.set_ylabels('Percent activated true-3d')
# %% calculate FASR complex
imdata.reset_index(inplace=True, drop=True)


def recruitment_cost(data, activated=1):
    """Calculate recruitment cost for each inner. :param activated: proportion
    of on-target fibers activated.

    Recruitment cost is defined as the ratio of number of stimulated
    off-target fibers to total number of off-target fibers. From
    https://iopscience.iop.org/article/10.1088/1741-2560/10/3/036010
    """
    for inner in pd.unique(data["inner"]):
        # get threshold for inner
        inner_data = data.query(f"inner == {inner}")["threshold"]
        assert len(data) > 200 & len(data) < 250
        # assert len(data['inner'])==len(set(data['inner']))
        # inner_thresh = np.amax(inner_data)
        # above line assumes 100% activation of on-target fibers, instead use activated
        inner_thresh = np.percentile(inner_data, activated * 100, method="higher")
        # get all off-target fiber thresholds
        off_thresh = data.query(f"inner != '{inner}'")["threshold"]
        # calculate recruitment cost
        cost = np.sum(off_thresh <= inner_thresh) / len(off_thresh)
        data.loc[data["inner"] == inner, "RC"] = cost
        # fascicle selectivity ratio is 1-RC
        data.loc[data["inner"] == inner, "FASR"] = 1 - cost
        fasr_dict = {
            "active_src_index": data["active_src_index"].iloc[0],
            "fiber_diam": data["fiber_diam"].iloc[0],
            "type": data["type"].iloc[0],
            "inner": inner,
            "RC": cost,
            "FASR": 1 - cost,
            "nerve_label": data["nerve_label"].iloc[0],
        }
        yield fasr_dict


imdatfasr = []
for contact_config in pd.unique(imdata["active_src_index"]):
    for fiber_diam in pd.unique(imdata["fiber_diam"]):
        for t in pd.unique(imdata["type"]):
            for nerve_label in pd.unique(imdata["nerve_label"]):
                imdatain = imdata.query(
                    f'active_src_index == "{contact_config}" and fiber_diam == {fiber_diam} and type == "{t}" and nerve_label=="{nerve_label}"'
                )
                imdatfasr.extend(recruitment_cost(imdatain, activated=0.90))

imdatfasr = pd.DataFrame(imdatfasr)
# %% plot FASR
fasrdiam = 3
sns.set(font_scale=2, style="whitegrid")
imdatfnewnonmerge = imdatfasr.query("fiber_diam in [@fasrdiam]")
imdatfnew = datamatch_merge(
    imdatfnewnonmerge.query('type=="extrusion"'),
    imdatfnewnonmerge.query('type=="true-3D"'),
    "RC",
    merge_cols=["active_src_index", "inner", "nerve_label"],
).drop(columns="type")
imdatfnew["RC-diff"] = imdatfnew["RC3d"] - imdatfnew["RC"]
for nerve_label in pd.unique(imdatfasr["nerve_label"]):
    plt.figure()
    imdatplot = imdatfnewnonmerge.query(f'nerve_label=="{nerve_label}"')
    imdatplot.RC *= 100
    # g = sns.boxplot(
    #     data=imdatplot.query('fiber_diam in [3]'),
    #     x='active_src_index',
    #     y='FASR-diff',
    #     # sharey=False,
    #     # hue='inner',
    #     # palette='rainbow',
    #     boxprops=dict(facecolor='none'),
    #     linewidth=3,
    #     # legend=False,
    #     whis=(0,100),
    # )
    g = sns.stripplot(
        data=imdatplot.query('fiber_diam in [@fasrdiam] and type=="extrusion"'),
        x="active_src_index",
        y="RC",
        hue="inner",
        palette="rainbow",
        legend=False,
        marker="s",
        dodge=True,
        edgecolor="k",
        linewidth=1,
        s=6,
    )
    g = sns.stripplot(
        data=imdatplot.query('fiber_diam in [@fasrdiam] and type=="true-3D"'),
        x="active_src_index",
        y="RC",
        hue="inner",
        palette="rainbow",
        legend=False,
        marker="o",
        dodge=True,
        edgecolor="black",
        linewidth=1,
        s=6,
    )
    dodges = {'2L': 0.68, '3R': 0.72, '5R': 0.73, '6R': 0.69}
    sns.pointplot(
        data=imdatplot.query("fiber_diam in [@fasrdiam]"),
        x="active_src_index",
        y="RC",
        hue="inner",
        palette="rainbow",
        legend=False,
        marker=None,
        err_kws={"linewidth": 2},
        dodge=dodges[nerve_label],
        linestyle="none",
        errorbar=("pi", 100),
    )
    # vertical line dashed in between each x value
    for i in range(5 if addln else 4):
        plt.axvline(i + 0.5, color="black", ls="--")
    # add legend for the two types. Create handles manually (gray marker with black outline)
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            label="extrusion",
            markerfacecolor="gray",
            markersize=10,
            markeredgewidth=1,
            markeredgecolor="black",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="true-3D",
            markerfacecolor="gray",
            markersize=10,
            markeredgewidth=1,
            markeredgecolor="black",
        ),
    ]
    legend_labels = ["extrusion", "true-3D"]
    # place legend above plot in 2 columns

    plt.legend(handles=legend_elements, labels=legend_labels, loc=(-0.1, 1.15), ncol=2)
    # plt.axhline(0,color='black',ls='--')
    plt.ylabel("off-target activated (%)")
    plt.ylim(0, 100)
    plt.xlabel("Active contact")
    plt.xticks(range(5), list(range(4)) + ["LN"])
    plt.title(f"Nerve: {nerve_label} - D: {fasrdiam} μm")
    plt.figure()
    # also plot ind
    # sys.exit()
# TODO: do FASR for line through nerve separation as well as for less than 100% activations
# %% MCT selectivity just the numbers
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(font_scale=1, style="white", context='paper')
import matplotlib.patheffects as PathEffects


def recruitment_cost_inner(data, activated=1, targetcol="inner"):
    """Calculate recruitment cost for each inner. :param activated:
    proportion of on-target fibers activated.

    Recruitment cost is defined as the ratio of number of stimulated
    off-target fibers to total number of off-target fibers. From
    https://iopscience.iop.org/article/10.1088/1741-2560/10/3/036010
    """
    for inner in pd.unique(data[targetcol]):
        # get threshold for inner
        inner_data = data.query(f"{targetcol} == {inner}")["threshold"]
        assert len(data) > 200 & len(data) < 250
        # assert len(data['inner'])==len(set(data['inner']))
        # inner_thresh = np.amax(inner_data)
        # above line assumes 100% activation of on-target fibers, instead use activated
        inner_thresh = np.percentile(inner_data, activated * 100, method="higher")
        # get all off-target fiber thresholds
        off_thresh = data.query(f"{targetcol} != '{inner}'")["threshold"]
        # calculate recruitment cost
        cost = np.sum(off_thresh <= inner_thresh) / len(off_thresh)
        data.loc[data[targetcol] == inner, "RC"] = cost
        # fascicle selectivity ratio is 1-RC
        data.loc[data[targetcol] == inner, "FASR"] = 1 - cost
        fasr_dict = {
            "active_src_index": data["active_src_index"].iloc[0],
            "fiber_diam": data["fiber_diam"].iloc[0],
            "type": data["type"].iloc[0],
            targetcol: inner,
            "RC": cost,
            "FASR": 1 - cost,
            "nerve_label": data["nerve_label"].iloc[0],
            "percent_ontarget": activated,
        }
        yield fasr_dict


lim = [-2000, 2000]
max_thk = 1000
allselectdata = []
for analysisdiam in imdata.fiber_diam.unique():
    print(analysisdiam)
    # sns.set(font_scale=1,style='white')
    for nerve_label, samplenum, r_cuff_in_pre_MCT in zip(
        ["2L", "3R", "5R", "6R"], [25212, 37212, 57212, 67212], [1000, 1500, 1000, 1000]
    ):
        imdata.reset_index(inplace=True, drop=True)
        imdatfasr = []
        for contact_config in pd.unique(imdata["active_src_index"]):
            for t in pd.unique(imdata["type"]):
                for percent_ontarget in [0.1, 0.5, 0.9]:
                    imdatain = imdata.query(
                        f'active_src_index == "{contact_config}" and fiber_diam == {analysisdiam} and type == "{t}" and nerve_label=="{nerve_label}"'
                    )
                    imdatain["percent_ontarget"] = percent_ontarget
                    imdatfasr.extend(recruitment_cost_inner(imdatain, activated=percent_ontarget, targetcol="inner"))
        imdatfasr = pd.DataFrame(imdatfasr)
        # now new plot, copy as above where extrusion and true-3D are plotted on the same graph with different markers, and the pointplot is added
        imdatplot = imdatfasr.query(f'fiber_diam in [{analysisdiam}]')
        imdatplot['RC'] *= 100
        imdatplot['fiber_diam'] = analysisdiam
        imdatplot['nerve_label'] = nerve_label
        allselectdata.append(imdatplot)
# %%
fiinalsel = datamatch_merge(
    pd.concat(allselectdata).query('type=="extrusion"'),
    pd.concat(allselectdata).query('type=="true-3D"'),
    "RC",
    merge_cols=["active_src_index", "fiber_diam", "inner", "nerve_label", "percent_ontarget"],
).drop(columns="type")
fiinalsel['resid'] = fiinalsel.RC3d - fiinalsel.RC
g = sns.FacetGrid(
    data=fiinalsel,
    row='fiber_diam',
)
g.map_dataframe(
    sns.violinplot,
    x='active_src_index',
    y='resid',
    hue='percent_ontarget',
    palette='plasma',
    legend=False,
    dodge=True,
    # jitter=False,
    linewidth=1,
    alpha=0.6,
    zorder=-2,
)
g.map_dataframe(
    sns.pointplot,
    x='active_src_index',
    y='resid',
    hue='percent_ontarget',
    palette='plasma',
    legend=False,
    estimator='median',
    errorbar=None,
    markeredgewidth=1,
    markeredgecolor='k',
    marker='s',
)

plt.gcf().set_size_inches(8, 8)
plt.ylim(-100, 100)
g.set_titles(row_template='D: {row_name} μm')
g.set_ylabels('Off target activation (%) \n(true-3D minus extrusion)')
g.set_xlabels('Active Contact')
for ax in g.axes.ravel():
    for pos, col in zip(np.arange(0.5, 4 if addln else 3, 1), ['gray', 'gray', 'gray', 'k']):
        ax.axvline(pos, ls='--', color=col)
# now absolute
fiinalsel['absresid'] = np.abs(fiinalsel.RC3d - fiinalsel.RC)
g = sns.FacetGrid(
    data=fiinalsel,
    row='fiber_diam',
)
g.map_dataframe(
    sns.violinplot,
    x='active_src_index',
    y='absresid',
    hue='percent_ontarget',
    palette='plasma',
    legend=False,
    dodge=True,
    # jitter=False,
    linewidth=1,
    alpha=0.6,
    zorder=-2,
)

g.map_dataframe(
    sns.pointplot,
    x='active_src_index',
    y='absresid',
    hue='percent_ontarget',
    palette='plasma',
    legend=False,
    estimator='median',
    errorbar=None,
    markeredgewidth=1,
    markeredgecolor='k',
    marker='s',
    dodge=0.53,
)

plt.gcf().set_size_inches(8, 4)
plt.ylim(0, 100)
g.set_titles(row_template='D: {row_name} μm')
g.set_ylabels('Off target activation (%) \n|true-3D minus extrusion|')
g.set_xlabels('Active Contact')
for ax in g.axes.ravel():
    for pos, col in zip(np.arange(0.5, 4, 1), ['gray', 'gray', 'gray', 'k']):
        ax.axvline(pos, ls='--', color=col)
simipledat = fiinalsel.copy()
simipledat['active_src_index'] = simipledat['active_src_index'].replace(
    {'1': 'MCT', '2': 'MCT', '3': 'MCT', '0': 'MCT'}
)
print(simipledat.groupby(['fiber_diam', 'active_src_index'])['resid', 'absresid'].median())
print(simipledat.groupby(['fiber_diam', 'active_src_index'])['resid', 'absresid'].min())
print(simipledat.groupby(['fiber_diam', 'active_src_index'])['resid', 'absresid'].max())

# lastly, calculate the min off target activated for each nerve, fiber_diam, active_src_index, and percent_ontarget
selectivedat = simipledat.groupby(['nerve_label', 'fiber_diam', 'active_src_index', 'percent_ontarget', 'inner'])[
    'RC3d'
].min()
# plot
g = sns.FacetGrid(
    data=selectivedat.reset_index(),
    row='fiber_diam',
)
g.map_dataframe(
    sns.barplot,
    x='active_src_index',
    y='RC3d',
    hue='percent_ontarget',
    palette='plasma',
    dodge=True,
    errorbar=None,
)
g.map_dataframe(
    sns.stripplot,
    x='active_src_index',
    color='k',
    y='RC3d',
    hue='percent_ontarget',
    dodge=True,
    marker='o',
    s=5,
    alpha=0.6,
    edgecolor='white',
    linewidth=1,
)
g.set_titles(row_template='D: {row_name} μm')
g.set_ylabels('Min. off target (%)')
g.set_xlabels('Cuff')
for ax in g.axes.ravel():
    ax.set_xticks([0, 1], ['LivaNova', 'MultiContact'])
plt.gcf().set_size_inches(2.5, 4)
# add values to barplots
# for ax in g.axes.ravel(): #uncomment to add numbers
#     for i in ax.containers:
#         ax.bar_label(i, label_type='edge',fmt='%.1f')
# %%
# plot
g = sns.FacetGrid(
    data=selectivedat.reset_index(),
    row='fiber_diam',
)
g.map_dataframe(
    sns.barplot,
    x='active_src_index',
    y='RC',
    hue='percent_ontarget',
    palette='plasma',
    dodge=True,
    errorbar=None,
)
g.map_dataframe(
    sns.stripplot,
    x='active_src_index',
    color='k',
    y='RC',
    hue='percent_ontarget',
    dodge=True,
    marker='o',
    s=5,
    alpha=0.6,
    edgecolor='white',
    linewidth=1,
)
g.set_titles(row_template='D: {row_name} μm')
g.set_ylabels('Min. off target (%)')
g.set_xlabels('Cuff')
for ax in g.axes.ravel():
    ax.set_xticks([0, 1], ['LivaNova', 'MultiContact'])
plt.gcf().set_size_inches(2.5, 4)
# add values to barplots
for ax in g.axes.ravel():
    for i in ax.containers:
        ax.bar_label(i, label_type='edge', fmt='%.1f')
# %%
# now absolute
fiinalsel['absresid'] = np.abs(fiinalsel.RC3d - fiinalsel.RC)
g = sns.FacetGrid(
    data=fiinalsel.query('active_src_index!="LN"'),
    row='fiber_diam',
)
g.map_dataframe(
    sns.violinplot,
    x='active_src_index',
    y='absresid',
    hue='percent_ontarget',
    palette='plasma',
    legend=False,
    dodge=True,
    # jitter=False,
    linewidth=1,
    alpha=0.6,
    zorder=-2,
    inner='quart',
)
g.map_dataframe(
    sns.swarmplot,
    x='active_src_index',
    y='absresid',
    hue='percent_ontarget',
    palette=['black'] * 3,
    legend=False,
    dodge=True,
    # jitter=0.2,
    s=2,
)
# g.map_dataframe(
#     sns.pointplot,
#     x='active_src_index',
#     y='absresid',
#     hue='percent_ontarget',
#     palette='plasma',
#     legend=False,
#     estimator='median',
#     errorbar=None,
#     markeredgewidth=1,
#     markeredgecolor='k',
#     marker='s',
#     dodge=0.53
# )
plt.gcf().set_size_inches(8, 4)
plt.ylim(0, 50)
g.set_titles(row_template='D: {row_name} μm')
g.set_ylabels('Off target activation (%) \n|true-3D minus extrusion|')
g.set_xlabels('Active Contact')
for ax in g.axes.ravel():
    for pos, col in zip(np.arange(0.5, 4, 1), ['gray', 'gray', 'gray', 'k']):
        ax.axvline(pos, ls='--', color=col)
simipledat = fiinalsel.copy()
simipledat['active_src_index'] = simipledat['active_src_index'].replace(
    {'1': 'MCT', '2': 'MCT', '3': 'MCT', '0': 'MCT'}
)
print(simipledat.groupby(['fiber_diam', 'active_src_index'])['resid', 'absresid'].median())
print(simipledat.groupby(['fiber_diam', 'active_src_index'])['resid', 'absresid'].min())
print(simipledat.groupby(['fiber_diam', 'active_src_index'])['resid', 'absresid'].max())
