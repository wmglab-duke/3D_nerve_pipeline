import json
import os

import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats
from scipy.stats import pearsonr, sem, variation

os.chdir('../../')
import sys

sys.path.pop(-2)
from src.core.plotter import datamatch_merge

mpl.rcParams['figure.dpi'] = 400
sns.set_style('whitegrid')
pd.options.mode.chained_assignment = None


def concordance_correlation_coefficient(y_true, y_pred):
    """Concordance correlation coefficient."""
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    # Raw data
    dct = {'y_true': y_true, 'y_pred': y_pred}
    df = pd.DataFrame(dct)
    # Remove NaNs
    df = df.dropna()
    # Pearson product-moment correlation coefficients
    y_true = df['y_true']
    y_pred = df['y_pred']
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
    for nsim in pd.unique(data['nsim']):
        # ax.set_title(f'fiber diam: {s}μm')
        corr = {}
        for sample in pd.unique(data['sample']):
            thisdat = data[(data["nsim"] == nsim) & (data["sample"] == sample)]
            corr[sample] = round(pearsonr(thisdat[comparison[0]], thisdat[comparison[1]])[0], 3)
        corrs.append(corr)
    return corrs


def addpwfd(data, sim, infile='plotconfig'):
    with open(f'examples/analysis/{infile}.json') as f:
        config = json.load(f)
    nsim_key = config['sim_data'][sim]['nsim_key']
    for nsim in nsim_key:
        pulse_width = nsim_key[nsim]['pulse_width']
        fiber_diam = nsim_key[nsim]['fiber_diam']
        nsim = int(nsim)
        data.loc[data.nsim == nsim, 'pulse_width'] = pulse_width
        data.loc[data.nsim == nsim, 'fiber_diam'] = fiber_diam
    return data


def pe(correct, est):
    """Calculate the percent error.

    :param correct: correct value
    :param est: estimated value
    :return: percent error
    """
    return 100 * abs(est - correct) / correct


def pe_noabs(correct, est):
    """Calculate the percent error.

    :param correct: correct value
    :param est: estimated value
    :return: percent error
    """
    return 100 * (est - correct) / correct


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

cath_comparison = ["cathodic", '3D', 3]
center_comparison = ["center", '3D', 3]
an_comparison = ["anodic", '3D', 3]
main_comparison = center_comparison

simNUM = input('Gimme simNUM: ')
if simNUM == 10:
    main_comparison = cath_comparison

pal2d3d = ['#d95f02', '#7570b3']

# code whole file to optionally run on 3 or 10 so that I can run all monpolar data
# %% set which comparison to run
print('threshload')
# base data
threshload = pd.read_csv(f'thresh_unmatched_sim{simNUM}_og.csv')

center = threshload['sample'].astype(str).str[2] == '1'
threshload.loc[center, 'contact'] = 'center'
threedcontact = threshload['sample'].astype(str).str[2] == '3'
threshload.loc[threedcontact, 'contact'] = '3D'
# set all inners and outers where 'type' is 3D to 0
threshload.loc[(threshload['type'] == '3D'), 'inner'] = 0
threshload.loc[(threshload['type'] == '3D'), 'outer'] = 0
# %% inners need to match the cathodic leading contact
# Query the rows with type '3D'
df_3d = threshload.query("type == '3D'")

# Query the rows with type '2D' and sample ending in '1'
df_2d = threshload.query(f"type == '2D' and contact in {main_comparison}")

# Merge the 3D and 2D data, keeping track of original row indices
merged_df = pd.merge(
    df_3d, df_2d, on=['nerve_label', 'master_fiber_index', 'nsim'], suffixes=('_3d', '_2d'), how='left'
)  # TODO remove this how

# Update the 'inner', 'outer', and 'fiber' columns in the original DataFrame
threshload.loc[df_3d.index, 'inner'] = merged_df['inner_2d'].values
threshload.loc[df_3d.index, 'outer'] = merged_df['outer_2d'].values
threshload.loc[df_3d.index, 'fiber'] = merged_df['fiber_2d'].values
# %%
if gogo == "initial":  # remove all where nerve_label length is >2
    threshload = threshload[threshload['nerve_label'].str.len() < 3]
threshload['type'] = threshload['type'].replace({'2D': 'extrusion', '3D': 'true-3D'})
threshload = addpwfd(threshload, str(simNUM), infile='plotconfig_og')
threshload['fiber_diam'] = threshload['fiber_diam'].astype(int)
# elif gogo=="deformedasc": #remove all where nerve_label does not contain "asc"
#     threshload = threshload[threshload['nerve_label'].str.contains("asc")]
#     #check the third digit of the sample number is 2, in that case, "contact" is cathodic. If 0, anodic
threshload['contact'] = (
    threshload['sample'].astype(str).str[2].replace({'0': 'anodic', '2': 'cathodic', '3': '3D', '1': 'center'})
)
threshload.sort_values(by='sample')

# %%Setup
print('matched')
# true-3D and extrusion trhresholds as different col same row
matched = datamatch_merge(
    threshload.query('type=="extrusion"'),
    threshload.query('type=="true-3D"'),
    'threshold',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='type')
# for matched boxplots
repeated = matched.melt(
    id_vars=[x for x in matched.columns if 'threshold' not in x],
    value_vars=['threshold3d', 'threshold'],
    var_name='type',
    value_name='threshold',
)
repeated['type'] = repeated['type'].replace({'threshold': 'extrusion', 'threshold3d': 'true-3D'})
repeated['type'] = pd.Categorical(repeated['type'], categories=['extrusion', 'true-3D'], ordered=True)
# %%
print('dose-response')
# dose response
drdat = threshload.copy().rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
drdat = calculate_dose_response(
    drdat,
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim'],
)
drdat.sort_values('modeltype', inplace=True)
drmatch = datamatch_merge(
    drdat.query('modeltype=="extrusion"'),
    drdat.query('modeltype=="true-3D"'),
    'percent_activated',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='modeltype')

# %% BEGIN DEFORMATION ANALYSIS
print('deformation')
# match fascicle thresholds
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshes = pd.concat(
    [pd.read_csv(f'thresh_unmatched_sim{simNUM}_og.csv'), pd.read_csv(f'thresh_unmatched_sim{simNUM}_def.csv')]
)
# remove all rows where nerve label contains "asc" and "sample" contains "3"
threshes = threshes[~((threshes['nerve_label'].str.contains('asc')) & (threshes['sample'].astype(str).str[2] == '3'))]
# duplicate all samples where nerve label contains "def" and "sample" contains "3"
defdat = threshes[(threshes['nerve_label'].str.contains('def')) & (threshes['sample'].astype(str).str[2] == '3')]
# replacee all "1" with "9" in sample column
defdat['sample'] = defdat['sample'].astype(str).str.replace('1', '9').astype(int)
# replace all "def" with "asc" in nerve label column
defdat['nerve_label'] = defdat['nerve_label'].str.replace('def', 'asc')
# combine with original data
threshes = pd.concat([threshes, defdat])
threshes['deformation'] = None
# where nerve label contains "def" deformation is "Structural"
threshes.loc[threshes['nerve_label'].str.contains('def'), 'deformation'] = 'Structural'
# where nerve label contains "asc" deformation is "ASCENT"
threshes.loc[threshes['nerve_label'].str.contains('asc'), 'deformation'] = 'ASCENT'
# else deformation is "none"
threshes.loc[threshes['deformation'].isna(), 'deformation'] = "Undeformed"
# strip "def" and "asc" from nerve labels
threshes['nerve_label'] = threshes['nerve_label'].str.replace('def', '').str.replace('asc', '')

newdefdat = threshes.copy()

# remove all nerve_label = '6L'
newdefdat = newdefdat[newdefdat['nerve_label'] != '6L']
newdefdat = newdefdat[newdefdat['nerve_label'] != '2R']

newdefdat['type'] = newdefdat['type'].replace({'2D': 'extrusion', '3D': 'true-3D'})
newdefdat = addpwfd(newdefdat, simNUM, infile='plotconfig_og')
newdefdat['fiber_diam'] = newdefdat['fiber_diam'].astype(int)
# contact is cathodic if the third digit of sample int is 2, anodic if 0
newdefdat['contact'] = (
    newdefdat['sample'].astype(str).str[2].replace({'2': 'cathodic', '1': 'center', '0': 'anodic', '3': '3D'})
)

# set deformation as ordered categorical
newdefdat['deformation'] = pd.Categorical(
    newdefdat['deformation'], categories=['Undeformed', 'ASCENT', 'Structural'], ordered=True
)
newdefdat['nerve_label'] = pd.Categorical(newdefdat['nerve_label'], categories=['2L', '3R', '5R', '6R'], ordered=True)
newdefdat['deformed'] = newdefdat['deformation'] != "Undeformed"
deftomatch = newdefdat.copy()

# remove unused colors from palette
defpal = [sns.color_palette('colorblind')[ind] for ind in [0, 2, 3, 5]]
defdefcomp = newdefdat.query('type=="true-3D" and deformation != "ASCENT"')
defdefcomp['deformed'] = defdefcomp['deformation'] != "Undeformed"
defsamples = ["2L", "3R", "5R", "6R"]
# %% new idea analysis
concats = []
for deftype in ["Undeformed", "Structural"]:
    thismatch = deftomatch.query('deformation==@deftype')
    matchednow = datamatch_merge(
        thismatch.query('type=="extrusion"'),
        thismatch.query('type=="true-3D"'),
        'threshold',
        merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
    ).drop(columns='type')
    matchednow['deformation'] = deftype
    concats.append(matchednow)
concats = pd.concat(concats)
concats['deformation'] = pd.Categorical(concats['deformation'], categories=['Undeformed', 'Structural'], ordered=True)

# %% dose-response
print('deformdr')
defdr = newdefdat.copy().rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
defdr = calculate_dose_response(
    defdr,
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim', 'deformation', 'pulse_width'],
)
sys.exit('prepdone')
# %%violinplot
newdef = newdefdat.copy()
newdef['deformation'] = pd.Categorical(
    newdef['deformation'], categories=['Undeformed', 'ASCENT', 'Structural'], ordered=True
)
plt.show()
g = sns.catplot(
    kind='violin',
    row='fiber_diam',
    data=newdef.query("fiber_diam in [3,13]"),
    # hue='type',
    x='type',
    y='threshold',
    sharey='row',
    errorbar='se',
    palette=pal2d3d,
    col='deformation',
    margin_titles=True,
)
# for ax in g.axes.ravel():
#     ax.set_yscale('log')
g.set_axis_labels('', 'Threshold (mA)')

# g.fig.set_size_inches(8,5)
# plt.yscale('log')
rownames(g, row_template="Threshold (mA)\nD: {row_name} μm")
g.set_titles(row_template='', col_template="{col_name}")
# %% dose-response
panelA = defdr.query(f"nerve_label in {defsamples} and contact in {main_comparison} and deformation=='Undeformed'")
panelA['panel'] = 'Undeformed extrusion\n(center slice) vs. true-3D'

panelB = defdr.query(f"modeltype=='true-3D' and deformation!='Structural' and contact in {main_comparison}")
panelB['panel'] = 'true-3D Undeformed\nvs. Deformed'

panelC = defdr.query(f"deformation!='Undeformed' and contact in {main_comparison}")
panelC['panel'] = 'Deformed extrusion\n(center slice) vs. true-3D'

panelD = defdr.query(f"deformation!='Undeformed' and contact in {cath_comparison}")
panelD['panel'] = 'Deformed extrusion (cathodic-\nleading slice) vs. true-3D'

alldr = pd.concat([panelA, panelB, panelC, panelD]).query('nerve_label=="2L"')
alldr['nerve_label'] = pd.Categorical(alldr['nerve_label'], categories=['2L'])

alldr['deformation'] = pd.Categorical(
    alldr['deformation'], categories=['Undeformed', 'ASCENT', 'Structural'], ordered=True
)

sns.set(font_scale=1.75, style='whitegrid')

plt.figure()
g = sns.relplot(
    kind='line',
    style='deformation',
    data=alldr.query(f"fiber_diam in [3]"),
    y='percent_activated',
    x='threshold',
    units='nerve_label',
    hue='modeltype',
    palette=pal2d3d,
    estimator=None,
    linewidth=4,
    col_wrap=2,
    facet_kws={'sharex': True, 'margin_titles': False},
    col='panel',
)

g.legend.set_title('')

for ax in g.axes.ravel():
    ax.set_ylim([0, 1])
    ax.set_xlim([0, None])
g.set_ylabels("Proportion fibers active")
g.set_titles(col_template='{col_name}')
# plt.subplots_adjust(hspace=0.25)

g.set_xlabels('Threshold (mA)')
# change the line width for the legend
for line, l in zip(g.legend.get_lines(), g.legend.get_texts()):
    line.set_linewidth(4.0)
    # if l.get_text() in ['deformation', 'modeltype']:
    #     l.set_text('')
g.fig.set_size_inches(12, 12)
for ax in g.axes.ravel():
    for loc in [0.1, 0.5, 0.9]:
        ax.axhline(loc, color='black', linestyle='-.', alpha=0.5, linewidth=2)
# %% dose-response
alldr = defdr.query(f"nerve_label in {defsamples}").query("nerve_label=='6R'")
alldr['deformed'] = alldr['deformation'] != "Undeformed"
allstore = []
for contact in alldr.contact.unique():
    if contact == "3D":
        continue
    thisdat = alldr.query('contact in [@contact, "3D"]')
    thisdat['contact'] = contact
    allstore.append(thisdat)
alldr = pd.concat(allstore)
sns.set(font_scale=1.75, style='whitegrid')

plt.figure()
g = sns.relplot(
    kind='line',
    style='deformation',
    data=alldr.query(f"fiber_diam in [3]"),
    y='percent_activated',
    x='threshold',
    units='nerve_label',
    hue='modeltype',
    palette=pal2d3d,
    estimator=None,
    linewidth=4,
    facet_kws={'sharex': True, 'margin_titles': True},
    col='deformed',
    row='contact',
)

g.legend.set_title('')

for ax in g.axes.ravel():
    ax.set_ylim([0, 1])
    ax.set_xlim([0, None])
g.set_ylabels("Proportion fibers active")
g.set_titles(col_template='{col_name}')
# plt.subplots_adjust(hspace=0.25)

g.set_xlabels('Threshold (mA)')
# change the line width for the legend
for line, l in zip(g.legend.get_lines(), g.legend.get_texts()):
    line.set_linewidth(4.0)
    # if l.get_text() in ['deformation', 'modeltype']:
    #     l.set_text('')
g.fig.set_size_inches(12, 12)
for ax in g.axes.ravel():
    for loc in [0.1, 0.5, 0.9]:
        ax.axhline(loc, color='black', linestyle='-.', alpha=0.5, linewidth=2)
g.set_titles(row_template='')
g.axes[0][0].set_title('Undeformed')
g.axes[0][1].set_title('Deformed')
rownames(g, row_template='{row_name}-slice\nProportion fibers active')
# %% threshold vals
for stringdat in ['Undeformed', 'Structural']:
    # deformstate = deformation!="Undeformed"
    thisdat = datamatch_merge(
        newdefdat.query(f'type=="extrusion" and deformation=="{stringdat}" and contact in {main_comparison}'),
        newdefdat.query(f'type=="true-3D" and deformation=="{stringdat}" and contact in {main_comparison}'),
        'threshold',
        merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
    ).drop(columns='type')
    for diam in [3, 13]:
        dadata = thisdat.query(f'fiber_diam==@diam and contact in {main_comparison}')
        perc = sum(dadata.threshold > dadata.threshold3d) / len(dadata.threshold)
        print(f'Percent higher thresh for {diam}:', round(perc, 3))
        res = np.abs(dadata.threshold3d - dadata.threshold)
        print(f'Mean abs residual for {diam}:', round(np.mean(res), 3), 'sem:', round(sem(res), 3))
        print(
            f'Max,min,mean residual for {diam}:',
            round(np.max(res), 3),
            round(np.min(res), 3),
            round(np.mean(res), 3),
            '+-',
            round(sem(res), 3),
        )
        mean_simthresh = np.mean(np.concatenate([dadata.threshold, dadata.threshold3d]))
        print(f'Mean sim thresh for {diam}:', round(mean_simthresh, 3))
        # max, min, mean, sem as a percentage of mean_simthresh
        print(
            f'Max,min,mean,sem as a percentage of mean_simthresh for {diam}:',
            round(np.max(res) / mean_simthresh, 3),
            round(np.min(res) / mean_simthresh, 3),
            round(np.mean(res) / mean_simthresh, 3),
            '+-',
            round(sem(res) / mean_simthresh, 3),
        )
        res = dadata.threshold3d - dadata.threshold
print(f'Difference between means for {diam}:', round(np.mean(res), 3), 'sem:', round(sem(res), 3))
# %% dose-response onset sat deformed
peses = []
pemeans = []
onsets_sats = {}
for stringdat in ['Undeformed', 'ASCENT', 'Structural']:
    subdat = newdefdat.query(f"deformation=='{stringdat}' and contact in {main_comparison}")
    sns.set(font_scale=1.75)
    sns.set_style('whitegrid')
    plt.figure()
    levels = {
        'onset': 10,
        'saturation': 90,
    }
    grouped = subdat.groupby(['sample', 'fiber_diam', 'type', 'sim', 'nerve_label', 'model', 'nsim', 'deformation'])
    analysis = grouped.agg(
        {
            'threshold': [
                lambda x: np.percentile(x, q=levels['onset']),
                np.median,
                lambda x: np.percentile(x, q=levels['saturation']),
            ]
        }
    )
    analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
    analysis.rename(columns={'threshold_<lambda_0>': 'onset', 'threshold_<lambda_1>': 'saturation'}, inplace=True)
    analysis = analysis.reset_index()
    # combine onset, saturation, and half into one column with identifier
    compiled_data = analysis.melt(
        id_vars=['sample', 'fiber_diam', 'sim', 'type', 'nerve_label', 'model', 'nsim'],
        value_vars=['onset', 'saturation'],
        var_name='level',
        value_name='threshold',
    )

    # set up facetgrid with nsim as row and level as columns
    compiled_data.reset_index(inplace=True)
    # set fiber_diam to category
    compiled_data.type = compiled_data.type.astype('category')
    # remove all rows where modulus sample with 0 is 0
    compiled_data = compiled_data.query('sample % 10 != 0')
    # add a units column with unique number for each combination of fiber_diam and level
    compiled_data['units'] = compiled_data.groupby(['fiber_diam', 'level', 'nerve_label']).ngroup()
    compiled_data['fiber_diam'] = compiled_data['fiber_diam'].astype(int)
    compiled_data.dropna(inplace=True)

    # calculate percent error for sample onset and saturation as well as population onset and saturation
    pes = []
    for level in ['onset', 'saturation']:
        for sam in compiled_data.nerve_label.unique():
            for fiber_diam in compiled_data.fiber_diam.unique():
                a = compiled_data.query(
                    f"nerve_label == '{sam}' and level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam}"
                )['threshold'].values
                b = compiled_data.query(
                    f"nerve_label == '{sam}' and level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}"
                )['threshold'].values
                assert len(a) == len(b) == 1
                pe_res = pe(a[0], b[0])
                pes.append({'level': level, 'nerve_label': sam, 'fiber_diam': fiber_diam, 'pe': pe_res})
    pes = pd.DataFrame(pes)
    pes["deformation"] = stringdat
    peses.append(pes)
    print('Max 3 um', stringdat, np.amax(pes.query('fiber_diam==3').pe))
    print('Max 13 um', stringdat, np.amax(pes.query('fiber_diam==13').pe))

    # now calculate percent error for population onset and saturation
    pemean = []
    for level in ['onset', 'saturation']:
        for fiber_diam in compiled_data.fiber_diam.unique():
            a = compiled_data.query(f"level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam}")[
                'threshold'
            ].values
            b = compiled_data.query(f"level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}")[
                'threshold'
            ].values
            pe_res = pe(np.median(a), np.median(b))
            pemean.append({'level': level, 'fiber_diam': fiber_diam, 'pe': pe_res})

    pemean = pd.DataFrame(pemean)
    pemean["deformation"] = stringdat
    pemeans.append(pemean)
    print('Max 3 um median', stringdat, np.amax(pemean.query('fiber_diam==3').pe))
    print('Max 13 um median', stringdat, np.amax(pemean.query('fiber_diam==13').pe))

    onsets_sats[stringdat] = compiled_data

sns.set(font_scale=2, style='white')
allpes = pd.concat(peses)
g = sns.relplot(
    data=allpes,
    kind='line',
    col='deformation',
    style='level',
    y='pe',
    x='fiber_diam',
    markers=True,
    color='black',
    units='nerve_label',
    estimator=None,
    linewidth=2,
    markersize=8,
    # hue='nerve_label'
)
# replace labels
g.set_xlabels('Fiber Diameter (μm)')
g.set_titles(col_template='Deformation: {col_name}')
g.set_ylabels("Percent Difference (%)")
plt.xticks([3, 5, 7, 9, 11, 13])
plt.ylim([0, None])
ylimpe = g.axes[0][0].get_ylim()
g.legend.set_title('')
sns.move_legend(g, 'lower center', bbox_to_anchor=(0.5, 0.6), ncol=1)
for line in g.legend.get_lines():
    line.set_linewidth(2)
for handle in g.legend.legend_handles:
    handle.set_markersize(10)
plt.ylim(0, 50)
allpescat = allpes.query('fiber_diam in [3,13]').groupby(['level', 'fiber_diam', 'deformation']).max()
pemeancat = pd.concat(pemeans).query('fiber_diam in [3,13]')
allpesmin = allpes.query('fiber_diam in [3,13]').groupby(['level', 'fiber_diam', 'deformation']).min()
# %% dose-response onset sat deformed
peses = []
pemeans = []
onsets_sats = {}
for stringdat in ['Undeformed', 'ASCENT', 'Structural']:
    subdat = newdefdat.query(f"deformation=='{stringdat}' and contact in {main_comparison}")
    sns.set(font_scale=1.75)
    sns.set_style('whitegrid')
    plt.figure()
    levels = {
        'onset': 10,
        'half': 50,
        'saturation': 90,
    }
    grouped = subdat.groupby(['sample', 'fiber_diam', 'type', 'sim', 'nerve_label', 'model', 'nsim', 'deformation'])
    analysis = grouped.agg(
        {
            'threshold': [
                lambda x: np.percentile(x, q=levels['onset']),
                lambda x: np.percentile(x, q=levels['half']),
                lambda x: np.percentile(x, q=levels['saturation']),
            ]
        }
    )
    analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
    analysis.rename(
        columns={'threshold_<lambda_0>': 'onset', 'threshold_<lambda_1>': 'half', 'threshold_<lambda_2>': 'saturation'},
        inplace=True,
    )
    analysis = analysis.reset_index()
    # combine onset, saturation, and half into one column with identifier
    compiled_data = analysis.melt(
        id_vars=['sample', 'fiber_diam', 'sim', 'type', 'nerve_label', 'model', 'nsim'],
        value_vars=['onset', 'half', 'saturation'],
        var_name='level',
        value_name='threshold',
    )

    # set up facetgrid with nsim as row and level as columns
    compiled_data.reset_index(inplace=True)
    # set fiber_diam to category
    compiled_data.type = compiled_data.type.astype('category')
    # add a units column with unique number for each combination of fiber_diam and level
    compiled_data['units'] = compiled_data.groupby(['fiber_diam', 'level', 'nerve_label']).ngroup()
    compiled_data['fiber_diam'] = compiled_data['fiber_diam'].astype(int)
    compiled_data.dropna(inplace=True)

    # calculate percent error for sample onset and saturation as well as population onset and saturation
    pes = []
    for level in ['onset', 'half', 'saturation']:
        for sam in compiled_data.nerve_label.unique():
            for fiber_diam in compiled_data.fiber_diam.unique():
                a = compiled_data.query(
                    f"nerve_label == '{sam}' and level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam}"
                )['threshold'].values
                b = compiled_data.query(
                    f"nerve_label == '{sam}' and level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}"
                )['threshold'].values
                assert len(a) == len(b) == 1
                pe_res = pe_noabs(a[0], b[0])
                pes.append({'level': level, 'nerve_label': sam, 'fiber_diam': fiber_diam, 'pe': pe_res})
    pes = pd.DataFrame(pes)
    pes["deformation"] = stringdat
    peses.append(pes)
    print('Max 3 um', stringdat, np.amax(pes.query('fiber_diam==3').pe))
    print('Max 13 um', stringdat, np.amax(pes.query('fiber_diam==13').pe))

    # now calculate percent error for population onset and saturation
    pemean = []
    for level in ['onset', 'half', 'saturation']:
        for fiber_diam in compiled_data.fiber_diam.unique():
            a = compiled_data.query(f"level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam}")[
                'threshold'
            ].values
            b = compiled_data.query(f"level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}")[
                'threshold'
            ].values
            pe_res = pe_noabs(np.median(a), np.median(b))
            pemean.append({'level': level, 'fiber_diam': fiber_diam, 'pe': pe_res})

    pemean = pd.DataFrame(pemean)
    pemean["deformation"] = stringdat
    pemeans.append(pemean)
    print('Max 3 um median', stringdat, np.amax(pemean.query('fiber_diam==3').pe))
    print('Max 13 um median', stringdat, np.amax(pemean.query('fiber_diam==13').pe))

    onsets_sats[stringdat] = compiled_data

sns.set(font_scale=2, style='white')
allpes = pd.concat(peses)
g = sns.relplot(
    data=allpes,
    kind='line',
    style='deformation',
    col='level',
    y='pe',
    x='fiber_diam',
    markers=True,
    color='black',
    linewidth=2,
    markersize=8,
    # hue='nerve_label',
    errorbar=('se', 1),
)
# replace labels
g.set_xlabels('Fiber Diameter (μm)')
g.set_titles(col_template='{col_name}')
g.set_ylabels("Percent Difference (%)")
plt.xticks([3, 5, 7, 9, 11, 13])
plt.ylim([0, None])
ylimpe = g.axes[0][0].get_ylim()
sns.move_legend(g, 'lower center', bbox_to_anchor=(0.45, 0.95), ncol=3)
for line in g.legend.get_lines():
    line.set_linewidth(2)
for handle in g.legend.legend_handles:
    handle.set_markersize(10)
g.legend
plt.ylim(-10, 50)
allpescat = allpes.query('fiber_diam in [3,13]').groupby(['level', 'fiber_diam', 'deformation']).max()
pemeancat = pd.concat(pemeans).query('fiber_diam in [3,13]')
allpesmin = allpes.query('fiber_diam in [3,13]').groupby(['level', 'fiber_diam', 'deformation']).min()
for ax in g.axes.ravel():
    ax.axhline(0, color='blue', linestyle='--', alpha=0.5)
# %% change in threshold after deformation
diffs = []
# calculate the change in onset and saturation for each sample (compared between undeformed and Structural deformed)
for nerve_label in pd.unique(onsets_sats["Undeformed"].nerve_label):
    for fiber_diam in pd.unique(onsets_sats["Undeformed"].fiber_diam):
        for level in ['onset', 'saturation']:
            undeformed = (
                onsets_sats["Undeformed"]
                .query(
                    f"type=='true-3D' and nerve_label == '{nerve_label}' and fiber_diam == {fiber_diam} and level == '{level}'"
                )['threshold']
                .values[0]
            )
            deformed = (
                onsets_sats["Structural"]
                .query(
                    f"type=='true-3D' and nerve_label == '{nerve_label}' and fiber_diam == {fiber_diam} and level == '{level}'"
                )['threshold']
                .values[0]
            )
            # print(f"{nerve_label} {fiber_diam} {level} {deformed-undeformed/undeformed}")
            diffs.append(
                {
                    'nerve_label': nerve_label,
                    'fiber_diam': fiber_diam,
                    'level': level,
                    'diff': 100 * (deformed - undeformed) / undeformed,
                }
            )
diffs = pd.DataFrame(diffs)
# print mean, std, min, max
print('Mean difference', np.mean(diffs['diff']))
print('Std difference', np.std(diffs['diff']))
print('Min difference', np.min(diffs['diff']))
print('Max difference', np.max(diffs['diff']))
# %% organization compare type with recruitment order
# remove all samples which dont end with a 2
sns.reset_orig()
sns.set(font_scale=1.25, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
# threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]
# # subtract 1 from sample if type is true-3D
# threshdat['sample'] = threshdat.apply(lambda x: x['sample'] - 1 if x.type == "true-3D" else x['sample'], axis=1)
data = threshload.query(f"fiber_diam in [3,13] and contact in {main_comparison} and nerve_label=='5R'")

sns.reset_orig()
sns.set(font_scale=1.25, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
# plt.suptitle('thresholds compared between extrusion (cathodic) and true-3D')
g = sns.FacetGrid(
    data,
    row="nerve_label",
    col='fiber_diam',
    sharey=False,
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x='type',
    y='threshold',
    units='master_fiber_index',
    estimator=None,
    palette='rainbow',
    hue='inner',
    linewidth=1,
    alpha=0.25,
)
g.map_dataframe(
    sns.stripplot,
    jitter=False,
    linewidth=1,
    x='type',
    y='threshold',
    palette='rainbow',
    hue='inner',
)

plt.subplots_adjust(top=0.9)
g.fig.set_size_inches(5, 5)
g.set_titles(col_template='D: {col_name} μm', row_template='')
g.set_xlabels('')

# handles, labs = plt.gca().get_legend_handles_labels()
# plt.legend(title='inner', handles=handles[14:], labels=labs[14:], bbox_to_anchor=[1, 2])

g.set_ylabels('Threshold (mA)')
# plt.suptitle('thresholds compared between extrusion (cathodic) and true-3D')

# %% organization compare nsim
# remove all samples which dont end with a 2
sns.reset_orig()
sns.set(font_scale=1.25, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')].query('nerve_label=="5R"')
data = threshdat.query('fiber_diam in [3,13] and sim == 3')
# apply minmax normalization within sample and nsim
data['Threshold (mA)'] = data.groupby(['sample', 'nsim', 'type'])['threshold'].transform(
    lambda x: (x - x.min()) / (x.max() - x.min())
)
data.query(f'contact in {main_comparison}', inplace=True)
# set nsim as categorical
data['Fiber Diameter'] = data['nsim'].astype('category')
# nsim 0 is three and nsim 5 is thirteen
data['Fiber Diameter'].cat.rename_categories({0: '3', 5: '13'}, inplace=True)
g = sns.FacetGrid(data, col='type', sharey=False, margin_titles=True)
g.map_dataframe(
    sns.lineplot,
    x='Fiber Diameter',
    y='Threshold (mA)',
    units='master_fiber_index',
    estimator=None,
    color='k',
    alpha=0.25,
    palette='rainbow',
    hue='inner',
)
g.map_dataframe(
    sns.stripplot,
    linewidth=1,
    x='Fiber Diameter',
    y='Threshold (mA)',
    jitter=False,
    palette='rainbow',
    hue='inner',
)
# plt.subplots_adjust(top=0.9)
# plt.suptitle('min-max normalized thresholds compared between 3 um and 13 um thresholds (cath extrusionEM)')
g.set_titles(col_template='{col_name}', row_template='')
g.set_ylabels('Normalized Threshold')
# plt.savefig('matchsim.png', dpi=400)
# g.axes.ravel()[0].set_ylabel('Threshold (mA)\nModel: extrusion')
# g.axes.ravel()[6].set_ylabel('Threshold (mA)\nModel: true-3D')
g.fig.set_size_inches(5, 5)
# %% organizatino compare nsim
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshdat = threshload.query(f'contact in {main_comparison}')

scores = []
for nerve in pd.unique(threshdat['nerve_label']):
    for n in [3, 13]:
        shortdat = threshdat.query(f'nerve_label=="{nerve}" and fiber_diam=={n}')
        data2d = shortdat.query('type=="extrusion"').sort_values('threshold').master_fiber_index
        data3d = shortdat.query('type=="true-3D"').sort_values('threshold').master_fiber_index
        rc = compute_reorder_cost(list(data2d), list(data3d))
        scores.append({'sample': nerve, 'fiber_diam': n, 'score2d3d': rc})
plt.figure()
scoredat = pd.DataFrame(scores)
scoredat['fiber_diam'] = pd.Categorical(scoredat['fiber_diam'].astype(int), categories=[3, 13], ordered=True)
ax = sns.boxplot(data=scoredat, y='score2d3d', x='fiber_diam', boxprops={'facecolor': 'white'}, whis=10)
ax = sns.stripplot(data=scoredat, y='score2d3d', x='fiber_diam', color='black', size=5, jitter=0.25)
ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('Activation Reordering')
plt.ylim(0, 0.5)
print(scoredat.groupby(['fiber_diam']).median())
# plt.gcf().set_size_inches([3, 5])
# %% organizatino compare type
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshdat = threshload.query(f'contact in {main_comparison}')

scores = []
for nerve in pd.unique(threshdat['nerve_label']):
    for model in ['extrusion', 'true-3D']:
        shortdat = threshdat.query(f'nerve_label=="{nerve}" and type=="{model}"')
        datasmol = shortdat.query('nsim==0').sort_values('threshold').master_fiber_index
        databeeg = shortdat.query('nsim==5').sort_values('threshold').master_fiber_index
        rc = compute_reorder_cost(list(datasmol), list(databeeg))
        scores.append({'sample': nerve, 'type': model, 'scoresmolbeeg': rc})
plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.boxplot(data=scoredat, y='scoresmolbeeg', x='type', boxprops={'facecolor': 'white'}, whis=10)
ax = sns.stripplot(data=scoredat, y='scoresmolbeeg', x='type', color='black', size=5, jitter=0.25)
plt.ylabel('Activation Reordering')
plt.xlabel('')
plt.ylim(0, 0.5)
print(scoredat.groupby(['type']).median())
# plt.gcf().set_size_inches([3, 5])
# %% activation order calc
newdefdr = defdr.copy()
newdefdr = newdefdr.query("'6R' in nerve_label and deformation=='Structural'").sort_values('modeltype')
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
newdefdr['percent_activated2d'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newdefdr.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        'modeltype == "true-3D" and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim and master_fiber_index == @row.master_fiber_index and deformation==@row.deformation'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newdefdr.loc[row.Index, 'percent_activated3d'] = val
newdefdr['contact'] = newdefdr['contact'].replace(
    {
        '3D': "true-3D",
        "cathodic": "cathodic-leading\nextrusion",
        "anodic": "anodic-leading\nextrusion",
        "center": "center\nextrusion",
    }
)
newdefdr['contact'] = pd.Categorical(
    newdefdr['contact'],
    categories=["true-3D", "cathodic-leading\nextrusion", "center\nextrusion", "anodic-leading\nextrusion"],
    ordered=True,
)
# %% plot activation order
sns.set(font_scale=1.25, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    row='fiber_diam',
    data=newdefdr.query("fiber_diam in [3]"),
    y='percent_activated',
    x='contact',
    units='nerve_label',
    palette='plasma',
    hue='percent_activated3d',
    # hue='inner',
    estimator=None,
    linewidth=0,
    facet_kws={'margin_titles': True},
    s=25,
)
plt.subplots_adjust(top=0.87)
# plt.suptitle(stringdat, x=0.37)
g.set_titles(row_template='', col_template='{col_name}')
g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion fibers activated')
g.set_xlabels('')
g.legend.remove()
norm = plt.Normalize(0, 1)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])
plt.axvline(0.5, linestyle='--', color='black', alpha=0.5)
# Remove the legend and add a colorbar
g.figure.colorbar(
    sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.8, label='Proportion true-3D\nfibers activated', pad=0.1
).ax.yaxis.set_ticks_position('left')
g.fig.set_size_inches(14, 4)
for i, con in enumerate(newdefdr.contact.sort_values().unique()[1:]):
    shortdat = newdefdr.query('fiber_diam==3 and contact==@con')
    data2d = shortdat.sort_values('percent_activated').master_fiber_index
    data3d = shortdat.sort_values('percent_activated3d').master_fiber_index
    rc = compute_reorder_cost(list(data2d), list(data3d))
    g.axes[0][0].text(i + 0.77, 1.065, round(rc, 3))
g.axes[0][0].text(-0.8, 1.07, 'Activation reordering: ')

# %% recruitment cost:
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
datahere = deftomatch.query('deformation in ["Structural","Undeformed"]')
scores = []
for comp in [main_comparison, cath_comparison, an_comparison]:
    for deformation in datahere.deformation.unique():
        threshdat = datahere.query(f'contact in {comp} and deformation==@deformation')
        for nerve in pd.unique(threshdat['nerve_label']):
            for n in [3, 13]:

                shortdat = threshdat.query(f'nerve_label=="{nerve}" and fiber_diam=={n}')
                data2d = shortdat.query('type=="extrusion"').sort_values('threshold').master_fiber_index
                data3d = shortdat.query('type=="true-3D"').sort_values('threshold').master_fiber_index
                rc = compute_reorder_cost(list(data2d), list(data3d))
                scores.append(
                    {'sample': nerve, 'fiber_diam': n, 'score2d3d': rc, 'deformation': deformation, 'slice': comp[0]}
                )
scoredat = pd.DataFrame(scores)
scoredat['slice'] = scoredat['slice'].replace(
    {"cathodic": "cathodic-\nleading", "anodic": "anodic-\nleading", "center": "center"}
)
scoredat['slice'] = pd.Categorical(
    scoredat['slice'], categories=["cathodic-\nleading", "center", "anodic-\nleading"], ordered=True
)
scoredat['fiber_diam'] = pd.Categorical(scoredat['fiber_diam'].astype(int), categories=[3, 13], ordered=True)
g = sns.FacetGrid(data=scoredat, col='slice', margin_titles=True, row="fiber_diam")
g.map_dataframe(
    sns.stripplot, y='score2d3d', x='deformation', hue='sample', palette=defpal, size=10, **dict(jitter=False)
)
g.map_dataframe(sns.lineplot, y='score2d3d', x='deformation', hue='sample', palette=defpal)
# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('Activation Reordering')
g.set_titles(col_template='{col_name}', row_template='D: {row_name} μm')
g.set_ylabels('Activation Reordering')
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
plt.ylim(0, 0.6)
print(scoredat.groupby(['deformation']).median())
plt.subplots_adjust(wspace=0)
plt.xlim(-1, 2)
g.fig.set_size_inches(6, 6)
for ax in g.axes.ravel():
    plt.sca(ax)
    plt.xticks(rotation=45)
# %% recruitment cost:
sns.set(font_scale=1)
sns.set_style('whitegrid')
datahere = deftomatch.query('deformation in ["Structural","Undeformed"]')
mathere = concats.copy()
scores = []
for comp in [main_comparison, cath_comparison, an_comparison]:
    for deformation in datahere.deformation.unique():
        threshdat = datahere.query(f'contact in {comp} and deformation==@deformation')
        for nerve in pd.unique(threshdat['nerve_label']):
            for n in [3, 13]:
                shortdat = threshdat.query(f'nerve_label=="{nerve}" and fiber_diam=={n}')
                data2d = shortdat.query('type=="extrusion"').sort_values('threshold').master_fiber_index
                data3d = shortdat.query('type=="true-3D"').sort_values('threshold').master_fiber_index
                rc = compute_reorder_cost(list(data2d), list(data3d))
                concathere = mathere.query(
                    f'nerve_label=="{nerve}" and fiber_diam=={n}and contact in {comp} and deformation==@deformation'
                )
                data2d = concathere.threshold
                data3d = concathere.threshold3d
                ccc = concordance_correlation_coefficient(list(data2d), list(data3d))
                scores.append(
                    {
                        'sample': nerve,
                        'fiber_diam': n,
                        'score2d3d': rc,
                        'deformation': deformation,
                        'slice': comp[0],
                        'CCC': ccc,
                    }
                )
nervescores = []
for comp in [main_comparison, cath_comparison, an_comparison]:
    for deformation in datahere.deformation.unique():
        threshdat = datahere.query(f'contact in {comp} and deformation==@deformation')
        for n in [3, 13]:
            shortdat = threshdat.query(f'fiber_diam=={n}')
            data2d = shortdat.query('type=="extrusion"').sort_values('threshold').master_fiber_index
            data3d = shortdat.query('type=="true-3D"').sort_values('threshold').master_fiber_index
            concathere = mathere.query(f'fiber_diam=={n}and contact in {comp} and deformation==@deformation')
            data2d = concathere.threshold
            data3d = concathere.threshold3d
            ccc = concordance_correlation_coefficient(list(data2d), list(data3d))
            nervescores.append(
                {'sample': nerve, 'fiber_diam': n, 'deformation': deformation, 'slice': comp[0], 'CCC': ccc}
            )
scoredat = pd.DataFrame(scores)
scoredat['slice'] = scoredat['slice'].replace(
    {"cathodic": "cathodic-\nleading", "anodic": "anodic-\nleading", "center": "center"}
)
scoredat['slice'] = pd.Categorical(
    scoredat['slice'], categories=["cathodic-\nleading", "center", "anodic-\nleading"], ordered=True
)
scoredat['fiber_diam'] = pd.Categorical(scoredat['fiber_diam'].astype(int), categories=[3, 13], ordered=True)
g = sns.FacetGrid(data=scoredat, margin_titles=True, col="fiber_diam")
g.map_dataframe(sns.scatterplot, y='score2d3d', x='slice', style='deformation', hue='sample', palette=defpal, size=10)
g.map_dataframe(sns.lineplot, y='score2d3d', x='slice', style='deformation', hue='sample', palette=defpal)
# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('Activation Reordering')
g.set_titles(col_template='D: {col_name} μm')
g.set_ylabels('Activation Reordering')
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
plt.ylim(0, 0.6)
print(scoredat.groupby(['deformation']).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
# g.fig.set_size_inches(8,4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel('extrusion slice')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[9:], labels=labs[9:], bbox_to_anchor=[1.2, 1])
# barplot
g = sns.FacetGrid(data=scoredat, margin_titles=True, col="fiber_diam")
g.map_dataframe(sns.barplot, y='score2d3d', x='slice', hue='deformation', errorbar='se', palette='binary')
g.map_dataframe(
    sns.stripplot,
    y='score2d3d',
    x='slice',
    hue='deformation',
    dodge=True,
    color='black',
    edgecolor='white',
    linewidth=0.5,
)

# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('Activation Reordering')
g.set_titles(col_template='D: {col_name} μm')
g.set_ylabels('Activation Reordering')
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
plt.ylim(0, 0.6)
print(scoredat.groupby(['deformation']).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
# g.fig.set_size_inches(8,4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel('extrusion slice')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[2:], labels=labs[2:], bbox_to_anchor=[0.5, 1.4], ncol=2)
# barplot
g = sns.FacetGrid(data=scoredat, margin_titles=True, col="fiber_diam")
g.map_dataframe(sns.barplot, y='CCC', x='slice', hue='deformation', errorbar='se', palette='binary')
g.map_dataframe(
    sns.stripplot, y='CCC', x='slice', hue='deformation', dodge=True, color='black', edgecolor='white', linewidth=0.5
)

# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('CCC')
g.set_titles(col_template='D: {col_name} μm')
g.set_ylabels('CCC')
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
# plt.ylim(0, 0.6)
print(scoredat.groupby(['deformation']).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
# g.fig.set_size_inches(8,4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel('extrusion slice')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[2:], labels=labs[2:], bbox_to_anchor=[0.5, 1.4], ncol=2)
# %%
g = sns.FacetGrid(data=pd.DataFrame(nervescores), margin_titles=True, col="fiber_diam")
g.map_dataframe(sns.barplot, y='CCC', x='slice', hue='deformation', errorbar='se', palette='binary')

# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('CCC')
g.set_titles(col_template='D: {col_name} μm')
g.set_ylabels('CCC')
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
# plt.ylim(0, 0.6)
print(scoredat.groupby(['deformation']).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
# g.fig.set_size_inches(8,4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel('extrusion slice')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[2:], labels=labs[2:], bbox_to_anchor=[0.5, 1.4], ncol=2)

# %% remake the above with individual calls of histplot
sns.set(font_scale=1.75, style='whitegrid')
newthreshz = threshload.copy()
newthreshz['activation_zpos'] = newthreshz['activation_zpos'] / 10000
fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
for nerve_label in pd.unique(newthreshz.nerve_label):
    for ax, modeltype in zip(axs, ["extrusion", "true-3D"]):
        g = sns.histplot(
            data=newthreshz.query(
                f"fiber_diam in [3,13] and contact in {main_comparison} and type=='{modeltype}'"
            ).rename(columns={'nerve_label': 'Sample'}),
            y='activation_zpos',
            # hue='fiber_diam',
            hue='Sample',
            # facet_kws={'sharex': False},
            # kind='kde',
            # palette=[sns.color_palette('binary')[5], sns.color_palette('binary')[2]],
            common_norm=False,
            # legend=False,
            # multiple="fill",
            element='poly',
            fill=False,
            bins=np.arange(1.5, 3.6, 0.1),
            ax=ax,
        )
# %% sample thresholds with unity line
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
# usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
usedata = matched
nsimdata = usedata.query(f'fiber_diam in [3,13] and contact in {main_comparison}')
g = sns.relplot(
    data=nsimdata.rename(columns={'nerve_label': 'Sample'}),
    kind='scatter',
    col='fiber_diam',
    x='threshold',
    y='threshold3d',
    # hue='Sample',
    color='white',
    # s=20,
    palette='colorblind',
    facet_kws={'sharex': False, 'sharey': False},
    edgecolor='black',
    linewidth=1,
    alpha=1,
)
# g.map(sns.regplot, 'threshold', 'threshold3d', scatter=False, color='k', label='linear fit')
# plot one one line out to max of true-3D

for diam, pos, ax in zip([3, 13], (0.2, 0.8), g.axes.ravel()):
    rdata = nsimdata.query(f'fiber_diam=={diam}')
    r, p = pearsonr(rdata.threshold, rdata.threshold3d)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} μm: {r ** 2:.2f}')
    ax.set_title(f'{diam} μm')
    ax.plot([0, rdata.threshold.max()], [0, rdata.threshold.max()], '--k', linewidth=2, label='unity line')
    ax.set_xlabel('extrusion Threshold (mA)')
    ax.set_yticks(ax.get_xticks())
    ax.set_aspect('equal', 'box')
    ax.set_xlim([0, None])
    ax.set_ylim([0, ax.get_xlim()[1]])
g.axes.ravel()[0].set_ylabel('true-3D Threshold (mA)')
g.axes.ravel()[1].set_ylabel('')
plt.legend(loc='lower right')
g.fig.set_size_inches([9.5, 5])
# %% Percent Error
sns.reset_orig()
sns.set_style('whitegrid')
sns.set(font_scale=1.5, style='whitegrid')
# apply pe to all rows of dataframe matched, with threshold3d as the correct value and threshold as the estimated value
matched['pe'] = matched.apply(lambda row: pe(row['threshold3d'], row['threshold']), axis=1)
plt.figure()
sns.barplot(data=matched, x='nerve_label', y='pe', hue='fiber_diam', errorbar='se', palette="RdPu")
# plt.title('Threshold Percent Error by sample and fiber diameter')
plt.legend(title='D (μm)', bbox_to_anchor=[0.55, 1.02], ncols=2)
plt.xlabel('')
plt.ylabel('Percent Difference (%)')
plt.gcf().set_size_inches([6, 5])

# calculate min, max, and mean percent error for each fiber diameter
pe_means = matched.groupby(['fiber_diam']).agg(np.mean)
pe_medians = matched.groupby(['fiber_diam']).agg(np.median)
pe_mins = matched.groupby(['fiber_diam']).agg(np.min)
pe_maxs = matched.groupby(['fiber_diam']).agg(np.max)
print("Percent Error by Fiber Diameter")
print("Mean: ", pe_means['pe'])
print("Median: ", pe_medians['pe'])
print("Min: ", pe_mins['pe'])
print("Max: ", pe_maxs['pe'])
# now do the same but for absolute error new column 'ae' is the absolute value of column 'pe'
matched['ae'] = matched['pe'].abs()
ae_means = matched.groupby(['fiber_diam']).agg(np.mean)
ae_medians = matched.groupby(['fiber_diam']).agg(np.median)
ae_mins = matched.groupby(['fiber_diam']).agg(np.min)
ae_maxs = matched.groupby(['fiber_diam']).agg(np.max)
print("Absolute Error by Fiber Diameter")
print("Mean: ", ae_means['ae'])
print("Median: ", ae_medians['pe'])
print("Min: ", ae_mins['ae'])
print("Max: ", ae_maxs['ae'])
# %% calculate mean, std, and median for each sample
outdata = (
    matched.query(f'contact in {main_comparison}')
    .groupby(['nerve_label', 'fiber_diam'])
    .agg({'threshold': ['mean', 'std', 'median'], 'threshold3d': ['mean', 'std', 'median']})
)
outdata.columns = ["_".join(col_name).rstrip('_') for col_name in outdata.columns]
outdata.reset_index(inplace=True)
outdata.dropna(inplace=True)
# calculate percent difference for mean, std, and median
outdata['mean_diff'] = (outdata.threshold_mean - outdata.threshold3d_mean) / outdata.threshold3d_mean * 100
outdata['std_diff'] = (outdata.threshold_std - outdata.threshold3d_std) / outdata.threshold3d_std * 100
outdata['median_diff'] = (outdata.threshold_median - outdata.threshold3d_median) / outdata.threshold3d_median * 100
# plot percent difference
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
sns.lineplot(data=outdata, x='fiber_diam', y='mean_diff', hue='nerve_label', palette='colorblind', marker='o')
plt.xticks(outdata.fiber_diam.unique())
plt.ylabel('Percent Difference')
plt.xlabel('Fiber Diameter (μm)')
plt.title('Mean Threshold Percent Difference')
plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0))
plt.show()
# %% threshold variances and coefficient of variation intrafascicle and inter
sns.set(font_scale=1.25, style='whitegrid')
vardat = repeated.query(f'contact in {main_comparison}')
grouped = vardat.groupby(['contact', 'sample', 'fiber_diam', 'type', 'inner'])
analysis = grouped.agg({'threshold': [np.var, np.mean, variation]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)

plt.figure()
g = sns.boxplot(
    data=analysis,
    y='threshold_variation',
    x='type',
    hue='fiber_diam',
    palette='RdPu',
)
plt.title('intrafascicle', pad=50)
plt.ylabel('Threshold Coefficient\nof Variation')
# plt.xlabel('Fiber Diameter (μm)')
# plt.yscale('log')
# plt.gca().get_legend().remove()
g.get_legend().set_title('D (μm)')
plt.gcf().set_size_inches([6, 4])
plt.title('')
plt.xlabel('')
# %% threshold variances and coefficient of variation intrafascicle and inter
sns.set(font_scale=1.75, style='whitegrid')
vardat = threshload.query(f'contact in {main_comparison}')
grouped = vardat.groupby(['contact', 'sample', 'fiber_diam', 'type', 'inner'])
analysis = grouped.agg({'peak_second_diff': [np.var, np.mean, variation]}).sort_values('type')
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)

plt.figure()
g = sns.boxplot(
    data=analysis,
    y='peak_second_diff_variation',
    x='type',
    hue='fiber_diam',
    palette='RdPu',
)
# plt.title('intrafascicle', pad=50)
plt.ylabel('Threshold CoV')
plt.xlabel('Fiber Diameter (μm)')
# plt.yscale('log')
plt.gca().get_legend().remove()
plt.gcf().set_size_inches([6, 4])
plt.yscale('log')
# %% imthera
imdata = pd.read_csv("thresh_unmatched_sim17_immy.csv")
imdata = addpwfd(imdata, '3')
imdata['contact'] = imdata['sample'].astype(str).str[2].replace({'2': 'cathodic', '0': 'anodic', '1': 'center'})
diams = {0: 3, 1: 13}
imdata['fiber_diam'] = imdata['fiberset_index'].replace(diams)

# inners need to match the center slice
innerid = {}
for inner in pd.unique(imdata.inner):
    innerid[inner] = pd.unique(imdata.query(f'inner == {inner} and type=="2D"')['master_fiber_index'])
    # get all rows where 3D and master fiber index is in the innerid
    imdata.loc[(imdata['type'] == '3D') & (imdata['master_fiber_index'].isin(innerid[inner])), 'inner'] = inner
imdata['type'] = imdata['type'].replace({'2D': 'extrusion', '3D': 'true-3D'})
im_matched = datamatch_merge(
    imdata.query('type=="extrusion"'),
    imdata.query('type=="true-3D"'),
    'threshold',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='type')
# %% imthera unity
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
# usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
g = sns.relplot(
    data=im_matched.rename(columns={'nerve_label': 'Sample'}).query('fiber_diam ==3'),
    kind='scatter',
    col='active_src_index',
    row='fiber_diam',
    x='threshold',
    y='threshold3d',
    # hue='Sample',
    color='white',
    # s=20,
    palette='rainbow',
    facet_kws={'sharex': False, 'sharey': 'row', 'margin_titles': True},
    edgecolor='black',
    linewidth=1,
    alpha=1,
    hue='inner',
    legend=False,
)

for ax in g.axes.ravel():
    lim = np.amax([ax.get_xlim()[1], ax.get_ylim()[1]])
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='unity line')
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])
    ax.set_aspect('equal')
rownames(g, row_template='true-3D Threshold (mA)\nD = {row_name} μm')
plt.legend(loc='lower right')
g.set_xlabels('extrusion Threshold (mA)')
g.set_titles(col_template='Active contact: {col_name}', row_template='')
plt.subplots_adjust(wspace=-0.6)

elceeceecee = []

for diam, contact in zip([3] * 6 + [13] * 6, [0, 1, 2, 3, 4, 5] * 2):
    rdata = im_matched.query(f'fiber_diam==@diam and active_src_index=={contact}')
    assert len(rdata) > 0
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    print(f'{contact} {diam} μm: {r ** 2:.2f}')
    elceeceecee.append(dict(diam=diam, contact=contact, value=r))
elcdata = pd.DataFrame(elceeceecee)
plt.figure()
sns.pointplot(data=elcdata, hue='diam', y='value', x='contact')
print(elcdata.groupby(['diam']).mean().round(2))
print(elcdata.groupby(['diam']).std().round(2))
# %% imthera contact config
sns.set(font_scale=1.5, style='whitegrid')
imdata['active_src_index'] = pd.Categorical(imdata['active_src_index'], categories=[0, 3, 1, 4, 2, 5], ordered=True)

g = sns.catplot(
    data=imdata.query('fiber_diam in [3]'),
    x='type',
    y='threshold',
    kind='strip',
    # row='fiber_diam',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    hue='inner',
    col='active_src_index',
    linewidth=0.75,
    s=30,
    dodge=True,
    col_wrap=2,
    palette='rainbow',
    edgecolor='k',
)
g.set_ylabels('Threshold (mA)')
g.set_xlabels('')
g.set_titles(col_template='Active contact: {col_name}', row_template='D: {row_name} μm')
g.fig.set_size_inches(10, 10)
# %% imthera contact config lineplot
sns.set(font_scale=1.5, style='whitegrid')
imdata['active_src_index'] = pd.Categorical(imdata['active_src_index'], categories=[0, 3, 1, 4, 2, 5], ordered=True)

g = sns.relplot(
    data=imdata.query('fiber_diam in [3]'),
    x='inner',
    y='threshold',
    kind='line',
    # row='fiber_diam',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    hue='type',
    col='active_src_index',
    linewidth=0.75,
    # s=30,
    # dodge=True,
    col_wrap=2,
    palette=pal2d3d,
    # edgecolor='k',
    errorbar=('pi', 100),
)
g.set_ylabels('Threshold (mA)')
g.set_xlabels('fascicle')
g.legend.set_title('')
g.set_titles(col_template='Active contact: {col_name}', row_template='D: {row_name} μm')
g.fig.set_size_inches(10, 10)
# %% activation order imthera
# dose response
imdr = imdata.copy()
imdr['active_src_index'] = pd.Categorical(imdata['active_src_index'], categories=[0, 1, 2, 3, 4, 5], ordered=True)

imdr = calculate_dose_response(
    imdr,
    'threshold',
    'percent_activated',
    grouping_columns=['type', 'sample', 'fiber_diam', 'sim', 'active_src_index'],
)
imdr.sort_values('type', inplace=True)

# now do
newim = imdr.copy()
newim = newim.sort_values('type')
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
newim['percent_activated2d'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newim.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newim.query(
        f'type == "true-3D" and fiber_diam == @row.fiber_diam and master_fiber_index == @row.master_fiber_index and active_src_index==@row.active_src_index'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newim.loc[row.Index, 'percent_activated3d'] = val
newim['type'] = pd.Categorical(newim['type'], ordered=True, categories=['true-3D', 'extrusion'])
# plot
sns.set(font_scale=1.5, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    data=newim.query('fiber_diam==3'),
    y='percent_activated',
    x='type',
    units='nerve_label',
    col='active_src_index',
    palette='plasma',
    hue='percent_activated3d',
    # hue='inner',
    estimator=None,
    linewidth=0,
    facet_kws={'margin_titles': True},
    s=25,
)
plt.subplots_adjust(top=0.87)
# plt.suptitle(stringdat, x=0.37)
g.set_titles(row_template='', col_template='Active contact: {col_name}')
g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion fibers activated')
g.set_xlabels('')
g.legend.remove()
norm = plt.Normalize(0, 1)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])

# Remove the legend and add a colorbar
g.figure.colorbar(
    sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.8, label='Proportion extrusion\nfibers activated', pad=0.025
).ax.yaxis.set_ticks_position('left')

g.fig.set_size_inches(18, 6)
# %% reorder cost for imthera
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

scores = []
for cc in range(6):
    shortdat = imdata.query(f'fiber_diam==3 and active_src_index=={cc}')
    data2d = shortdat.query('type=="extrusion"').sort_values('threshold').master_fiber_index
    data3d = shortdat.query('type=="true-3D"').sort_values('threshold').master_fiber_index
    rc = compute_reorder_cost(list(data2d), list(data3d))
    scores.append({'active_src_index': cc, 'score2d3d': rc})
plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.barplot(data=scoredat, y='score2d3d', x='active_src_index', color='k')
ax.set_xlabel('Active Contact')
plt.ylabel('Activation Reordering')
plt.ylim(0, 0.5)
# %% reorder cost for imthera
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

scores = []
for cc in range(6):
    shortdat = imdata.query(f'fiber_diam==3 and active_src_index=={cc}')
    data2d = shortdat.query('type=="extrusion"')['threshold']
    data3d = shortdat.query('type=="true-3D"')['threshold']
    ccc = concordance_correlation_coefficient(data3d, data2d)
    scores.append({'active_src_index': cc, 'score2d3d': ccc})
plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.barplot(data=scoredat, y='score2d3d', x='active_src_index', color='k')
ax.set_xlabel('Active Contact')
plt.ylabel('CCC')

# %% reorder cost for imthera fascicle specific
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

mean_im = imdata.groupby(['active_src_index', 'inner', 'type', 'fiber_diam']).mean().reset_index()

scores = []
for cc in range(6):
    shortdat = mean_im.query(f'fiber_diam==3 and active_src_index=={cc}')
    data2d = shortdat.query('type=="extrusion"').sort_values('threshold').inner
    data3d = shortdat.query('type=="true-3D"').sort_values('threshold').inner
    rc = compute_reorder_cost(list(data2d), list(data3d))
    scores.append({'active_src_index': cc, 'score2d3d': rc})
plt.figure()
ax = sns.barplot(data=scoredat, y='score2d3d', x='active_src_index', color='k')
ax.set_xlabel('Active Contact')
plt.ylabel('Activation Reordering')
plt.ylim(0, 0.5)

# %% dose-response imthera
sns.set_style('white')
imdr = imdata.copy().rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
imdr = calculate_dose_response(
    imdr,
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim', 'active_src_index'],
)
imdr.sort_values('modeltype', inplace=True)
imdrmatch = datamatch_merge(
    imdr.query('modeltype=="extrusion"'),
    imdr.query('modeltype=="true-3D"'),
    'percent_activated',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='modeltype')

fig, axs = plt.subplots(1, 6, sharey=True, sharex=True)
for acsrc, ax in zip(sorted(pd.unique(imdr.active_src_index)), axs):
    # plt.figure()
    sns.set(font_scale=1.5, style='white')
    g = sns.scatterplot(
        data=imdr.query(f"fiber_diam in [3] and active_src_index==@acsrc"),
        y='percent_activated',
        x='threshold',
        hue='inner',
        palette='rainbow',
        linewidth=0,
        ax=ax,
        # legend=False
    )
    sns.lineplot(
        data=imdr.query(f"fiber_diam in [3] and active_src_index==@acsrc"),
        y='percent_activated',
        x='threshold',
        color='k',
        linewidth=2,
        estimator=None,
        units='modeltype',
        legend=False,
        zorder=1,
        ax=ax,
        style='modeltype',
    )

    ax.set_xlabel('Threshold (mA)')
    ax.set_ylim([0, 1])
    # plt.gcf().set_size_inches(8,4)
    # create legend, circle = extrusion, X= true-3D
    # create handles manually
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], linestyle='-', color='k', label='extrusion', linewidth=2),
        Line2D([0], [0], linestyle='--', color='k', label='true-3D', linewidth=2),
    ]
    legend_labels = ['extrusion', 'true-3D']
    g.legend(handles=legend_elements, labels=legend_labels, loc='lower right')
    ax.set_title(f'Active contact: {acsrc}')
fig.set_size_inches(24, 4)
axs[0].set_ylabel('Proportion Activated')
plt.subplots_adjust(wspace=0)
# %% for each contact configuration, and fiber diameter, calculate the fascicle selectivity ratio for each inner
imdata.reset_index(inplace=True, drop=True)


def recruitment_cost(data):
    """Calculate recruitment cost for each inner.

    Recruitment cost is defined as the ratio of number of stimulated off-target fibers to total number of off-target fibers.
    From https://iopscience.iop.org/article/10.1088/1741-2560/10/3/036010
    """
    for inner in pd.unique(data["inner"]):
        # get threshold for inner
        inner_data = data.query(f"inner == {inner}")['threshold']
        assert len(data) > 200 & len(data) < 250
        inner_thresh = np.amax(inner_data)
        # get all off-target fiber thresholds
        off_thresh = data.query(f"inner != '{inner}'")['threshold']
        # calculate recruitment cost
        cost = np.sum(off_thresh <= inner_thresh) / len(off_thresh)
        data.loc[data['inner'] == inner, 'RC'] = cost
        # fascicle selectivity ratio is 1-RC
        data.loc[data['inner'] == inner, 'FASR'] = 1 - cost
        fasr_dict = {
            'active_src_index': data['active_src_index'].iloc[0],
            'fiber_diam': data['fiber_diam'].iloc[0],
            'type': data['type'].iloc[0],
            'inner': inner,
            'RC': cost,
            'FASR': 1 - cost,
        }
        yield fasr_dict


imdatfasr = []
for contact_config in pd.unique(imdata['active_src_index']):
    for fiber_diam in pd.unique(imdata['fiber_diam']):
        for t in pd.unique(imdata['type']):
            imdatain = imdata.query(
                f'active_src_index == {contact_config} and fiber_diam == {fiber_diam} and type == "{t}"'
            )
            imdatfasr.extend(recruitment_cost(imdatain))

imdatfasr = pd.DataFrame(imdatfasr)
# %% plot FASR
sns.set(font_scale=2, style='white')
g = sns.relplot(
    data=imdatfasr.query('fiber_diam in [3]'),
    x='active_src_index',
    y='FASR',
    kind='line',
    col='inner',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    hue='inner',
    style='type',
    palette='rainbow',
    linewidth=3,
    col_wrap=7,
)
plt.ylim(0, 1)
g.set_xlabels('Active contact')
g.set_titles(col_template='fascicle: {col_name}')
plt.xticks(range(7))
# g.fig.set_size_inches([8, 10])
g.legend.remove()
g.fig.legend(handles=g.legend.legendHandles[-2:], labels=['extrusion', 'true-3D'], bbox_to_anchor=[0.94, 0.23])
# %% plot FASR
sns.set(font_scale=2, style='whitegrid')
imdatfnew = imdatfasr.query('fiber_diam in [3]')
imdatfnew = datamatch_merge(
    imdatfnew.query('type=="extrusion"'),
    imdatfnew.query('type=="true-3D"'),
    'FASR',
    merge_cols=["active_src_index", "inner"],
).drop(columns='type')
imdatfnew['FASR-diff'] = imdatfnew['FASR3d'] - imdatfnew['FASR']
g = sns.boxplot(
    data=imdatfnew,
    x='active_src_index',
    y='FASR-diff',
    # sharey=False,
    # hue='inner',
    # palette='rainbow',
    boxprops=dict(facecolor='none'),
    linewidth=3,
    # legend=False,
    whis=(0, 100),
)
g = sns.lineplot(
    data=imdatfnew,
    x='active_src_index',
    y='FASR-diff',
    # sharey=False,
    hue='inner',
    palette='rainbow',
    linewidth=2,
    legend=False,
    alpha=0.5,
)
plt.ylim(-1, 1)
plt.xlabel('Active contact')
plt.xticks(range(6))
# %% plot FASR
sns.set(font_scale=1.75, style='white')
g = sns.relplot(
    data=imdatfasr.query('fiber_diam in [3]'),
    x='inner',
    y='FASR',
    kind='line',
    col='active_src_index',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    style='type',
    palette='rainbow',
    linewidth=3,
    col_wrap=3,
)
plt.ylim(0, 1)
# %% remake the above with individual calls of histplot
sns.set(font_scale=1.75, style='whitegrid')
newthreshz = newdefdat.query("deformation=='Structural'")
newthreshz['activation_zpos'] = newthreshz['activation_zpos'] / 10000
fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
for nerve_label in pd.unique(newthreshz.nerve_label):
    for ax, modeltype in zip(axs, ["extrusion", "true-3D"]):
        g = sns.histplot(
            data=newthreshz.query(
                f"fiber_diam in [3,13] and contact in {main_comparison} and nerve_label=='{nerve_label}' and type=='{modeltype}'"
            ).rename(columns={'nerve_label': 'Sample'}),
            y='activation_zpos',
            hue='fiber_diam',
            # hue='Sample',
            # facet_kws={'sharex': False},
            # kind='kde',
            palette=[sns.color_palette('binary')[5], sns.color_palette('binary')[2]],
            common_norm=False,
            # legend=False,
            # multiple="fill",
            element='poly',
            fill=False,
            bins=np.arange(1.5, 3.6, 0.1),
            ax=ax,
        )
# delete both legends and remake my own
for ax in axs:
    ax.get_legend().remove()
# make my own legend
axs[0].plot([], [], color=sns.color_palette('binary')[5], label='3 μm', linewidth=2)
axs[0].plot([], [], color=sns.color_palette('binary')[2], label='13 μm', linewidth=2)
# put legenbd to right of figure
axs[0].axhspan(2.005, 2.205, color='red', alpha=0.2, label='contacts')
axs[1].axhspan(2.005, 2.205, color='red', alpha=0.2, label='_')
axs[0].axhspan(2.805, 3.005, color='red', alpha=0.2, label='_')
axs[1].axhspan(2.805, 3.005, color='red', alpha=0.2, label='_')

axs[0].legend(loc='center left', bbox_to_anchor=(2.2, 0.5))
axs[0].set_xlim(reversed(axs[0].get_xlim()))
axs[0].set_ylabel("Activation Location (cm)")
axs[0].set_title("extrusion-100%")
axs[1].set_title("true-3D-100%")

# %% Correlationextrusiontrue-3D
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = matched.rename(columns={'threshold': 'thresholdextrusion'})
usedata.threshold3d = usedata.threshold3d.astype(float)
for comparison in [
    # ['thresholdextrusion', 'threshold3d'],
    # ['thresholdextrusion', 'peri_thk'],
    # ['threshold3d', 'peri_thk'],
    # ['thresholdextrusion', 'minimum_efib_distance'],
    ['thresholdextrusion', 'peak_second_diff'],
    # ['thresholdextrusion', 'peak_abs_second_diff'],
    # ['peak_second_diff', 'peri_thk'],
    # ['peak_second_diff', 'minimum_efib_distance'],
    # ['peak_second_z', 'activation_zpos'],
    # ['peak_second_diff_node', 'apnode'],
]:
    corrs = usedata.groupby(['sample', 'fiber_diam', 'contact', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
    corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    corrs['contact'] = pd.Categorical(corrs['contact'], categories=['cathodic', 'center', 'anodic'], ordered=True)
    plt.figure()
    # g = sns.FacetGrid(data=corrs,col='contact')
    # g.map_dataframe(sns.stripplot, hue='Sample', y='correlation',dodge=True, palette='colorblind')
    # g.map_dataframe(sns.boxplot,y='correlation', boxprops={'facecolor':'None'},whis=100)
    sns.swarmplot(data=corrs, palette='RdPu', y='contact', hue='fiber_diam', x='correlation', dodge=True, s=6)
    sns.boxplot(data=corrs, y='contact', x='correlation', boxprops={'facecolor': 'None'}, whis=100)
    # plt.subplots_adjust(top=0.8)
    plt.title(f'Correlation between \n{comparison[0]} and {comparison[1]}', pad=25)
    plt.gca().set_xlim([-1, 1])
    plt.legend(title="D (μm)", bbox_to_anchor=(1, 1))
    plt.ylabel('')
    plt.gca().set_yticklabels('')
    plt.gcf().set_size_inches([6, 5])
    # plt.yticks([0,1,2],['anodic-leading','center','cathodic-leading'])
    # plt.ylabel('extrusion slice')
    # plt.title('')

    # make lmplot to accompany
    sns.lmplot(
        data=usedata.query('fiber_diam in [3,13]'),
        x=comparison[0],
        y=comparison[1],
        hue='nerve_label',
        col='contact',
        col_order=['cathodic', 'center', 'anodic'],
        palette='colorblind',
        row='fiber_diam',
        row_order=[3, 13],
        facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
        scatter_kws={'linewidth': 1, 'edgecolor': 'k'},
    )

# %% Correlation 3Donly
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = threshload.query('contact=="3D"')
for comparison in [
    # ['threshold', 'tortuosity'],
    # ['threshold', 'nodal_tortuosity'],
    # ['threshold', 'cuff_tortuosity'],
    # ['threshold', 'peri_thk_act_site'],
    # ['threshold', 'smallest_thk_under_cuff'],
    # ['threshold','peri_thk_act_site'],
    # ['smallest_thk_under_cuff', 'peri_thk_act_site'],
    # # TODO try smallest thk over whole cuffspan as well as over both cuff spans
    # ['threshold', 'minimum_efib_distance'],
    ['threshold', 'peak_second_diff'],
    ['threshold', 'peak_abs_second_diff'],
    # ['peak_second_diff', 'peri_thk_act_site'],
    # ['peak_second_diff', 'tortuosity'],
    # ['peak_second_diff', 'nodal_tortuosity'],
    ['peak_second_z', 'activation_zpos'],
    # ['peak_second_long_pos', 'long_ap_pos'],
    # ['peak_second_diff_node', 'apnode'],
    # ['peak_second_diff', 'minimum_efib_distance'],
]:
    corrs = usedata.groupby(['sample', 'fiber_diam', 'contact', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
    corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    corrs['contact'] = pd.Categorical(corrs['contact'], ordered=True)
    plt.figure()
    # g = sns.FacetGrid(data=corrs,col='contact')
    # g.map_dataframe(sns.stripplot, hue='Sample', y='correlation',dodge=True, palette='colorblind')
    # g.map_dataframe(sns.boxplot,y='correlation', boxprops={'facecolor':'None'},whis=100)
    sns.swarmplot(data=corrs, palette='RdPu', y='contact', hue='fiber_diam', x='correlation', dodge=True, s=6)
    sns.boxplot(data=corrs, y='contact', x='correlation', boxprops={'facecolor': 'None'}, whis=100)
    # plt.subplots_adjust(top=0.8)
    plt.title(f'Correlation between \n{comparison[0]} and {comparison[1]}', pad=25)
    plt.gca().set_xlim([-1, 1])
    plt.legend(title="D (μm)", bbox_to_anchor=(1, 1))
    plt.ylabel('')
    plt.gca().set_yticklabels('')
    plt.gcf().set_size_inches([6, 5])
    # make lmplot to accompany
    sns.lmplot(
        data=usedata.query('fiber_diam in [3,13]'),
        x=comparison[0],
        y=comparison[1],
        hue='nerve_label',
        row='fiber_diam',
        row_order=[3, 13],
        palette='colorblind',
        facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
        scatter_kws={'linewidth': 1, 'edgecolor': 'k'},
    )
# %% Correlationextrusiontrue-3D deformed
threeddefmatch = deftomatch.query('deformation=="Structural"')
deffinalmatch = datamatch_merge(
    threeddefmatch.query('type=="extrusion"'),
    threeddefmatch.query('type=="true-3D"'),
    'threshold',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='type')
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = deffinalmatch.rename(columns={'threshold': 'thresholdextrusion'})
for comparison in [
    ['thresholdextrusion', 'threshold3d'],
    # ['thresholdextrusion', 'peri_thk'],
    # ['threshold3d', 'peri_thk'],
    # # ['thresholdextrusion', 'minimum_efib_distance'],
    # ['thresholdextrusion', 'peak_second_diff'],
    # ['peak_second_diff', 'peri_thk'],
    # # ['peak_second_diff', 'minimum_efib_distance'],
    # ['peak_second_z', 'activation_zpos'],
    # ['peak_second_diff_node', 'apnode'],
]:
    corrs = usedata.groupby(['sample', 'fiber_diam', 'contact', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
    corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    corrs['contact'] = pd.Categorical(corrs['contact'], categories=['cathodic', 'center', 'anodic'], ordered=True)
    plt.figure()
    # g = sns.FacetGrid(data=corrs,col='contact')
    # g.map_dataframe(sns.stripplot, hue='Sample', y='correlation',dodge=True, palette='colorblind')
    # g.map_dataframe(sns.boxplot,y='correlation', boxprops={'facecolor':'None'},whis=100)
    sns.swarmplot(data=corrs, palette='RdPu', y='contact', hue='fiber_diam', x='correlation', dodge=True, s=6)
    sns.boxplot(data=corrs, y='contact', x='correlation', boxprops={'facecolor': 'None'}, whis=100)
    # plt.subplots_adjust(top=0.8)
    plt.title(f'Correlation between \n{comparison[0]} and {comparison[1]}', pad=25)
    plt.gca().set_xlim([-1, 1])
    plt.legend(title="D (μm)", bbox_to_anchor=(1, 1))
    plt.ylabel('')
    plt.gca().set_yticklabels('')
    plt.gcf().set_size_inches([6, 5])

    # make lmplot to accompany
    sns.lmplot(
        data=usedata.query('fiber_diam in [3,13]'),
        x=comparison[0],
        y=comparison[1],
        hue='nerve_label',
        col='contact',
        col_order=['cathodic', 'center', 'anodic'],
        palette='colorblind',
        row='fiber_diam',
        row_order=[3, 13],
        facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
        scatter_kws={'linewidth': 1, 'edgecolor': 'k'},
    )
# %% Correlation 3Donly deformed
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = newdefdat.query('contact=="3D" and deformation=="Structural"')
for comparison in [
    # ['threshold', 'tortuosity'],
    # ['threshold', 'nodal_tortuosity'],
    ['threshold', 'cuff_tortuosity'],
    # ['threshold', 'peri_thk_act_site'],
    ['threshold', 'smallest_thk_under_cuff'],
    # ['threshold', 'peri_thk_act_site'],
    ['smallest_thk_under_cuff', 'peri_thk_act_site'],
    # TODO try smallest thk over whole cuffspan as well as over both cuff spans
    # ['threshold', 'minimum_efib_distance'],
    # ['threshold', 'peak_second_diff'],
    # ['peak_second_diff', 'peri_thk_act_site'],
    # ['peak_second_diff', 'tortuosity'],
    # ['peak_second_diff', 'nodal_tortuosity'],
    # ['peak_second_z', 'activation_zpos'],
    # ['peak_second_long_pos', 'long_ap_pos'],
    # ['peak_second_diff_node', 'apnode'],
    # ['peak_second_diff', 'minimum_efib_distance'],
]:
    corrs = usedata.groupby(['sample', 'fiber_diam', 'contact', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
    corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    corrs['contact'] = pd.Categorical(corrs['contact'], ordered=True)
    plt.figure()
    # g = sns.FacetGrid(data=corrs,col='contact')
    # g.map_dataframe(sns.stripplot, hue='Sample', y='correlation',dodge=True, palette='colorblind')
    # g.map_dataframe(sns.boxplot,y='correlation', boxprops={'facecolor':'None'},whis=100)
    sns.swarmplot(data=corrs, palette='RdPu', y='contact', hue='fiber_diam', x='correlation', dodge=True, s=6)
    sns.boxplot(data=corrs, y='contact', x='correlation', boxprops={'facecolor': 'None'}, whis=100)
    # plt.subplots_adjust(top=0.8)
    plt.title(f'Correlation between \n{comparison[0]} and {comparison[1]}', pad=25)
    plt.gca().set_xlim([-1, 1])
    plt.legend(title="D (μm)", bbox_to_anchor=(1, 1))
    plt.ylabel('')
    plt.gca().set_yticklabels('')
    plt.gcf().set_size_inches([6, 5])
    # make lmplot to accompany
    sns.lmplot(
        data=usedata.query('fiber_diam in [3,13]'),
        x=comparison[0],
        y=comparison[1],
        hue='nerve_label',
        row='fiber_diam',
        row_order=[3, 13],
        palette='colorblind',
        facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
        scatter_kws={'linewidth': 1, 'edgecolor': 'k'},
    )
# %% Percent Error deformed
sns.reset_orig()
sns.set_style('whitegrid')
sns.set(font_scale=1.5, style='whitegrid')
# apply pe to all rows of dataframe matched, with threshold3d as the correct value and threshold as the estimated value
deffinalmatch['pe'] = deffinalmatch.apply(lambda row: pe(row['threshold3d'], row['threshold']), axis=1)
plt.figure()
sns.barplot(data=deffinalmatch, x='nerve_label', y='pe', hue='fiber_diam', errorbar='se', palette="RdPu")
# plt.title('Threshold Percent Error by sample and fiber diameter')
plt.legend(title='D (μm)', bbox_to_anchor=[0.55, 1.02], ncols=2)
plt.xlabel('')
plt.ylabel('Percent Difference (%)')
plt.gcf().set_size_inches([6, 5])

# calculate min, max, and mean percent error for each fiber diameter
pe_means = deffinalmatch.groupby(['fiber_diam']).agg(np.mean)
pe_medians = deffinalmatch.groupby(['fiber_diam']).agg(np.median)
pe_mins = deffinalmatch.groupby(['fiber_diam']).agg(np.min)
pe_maxs = deffinalmatch.groupby(['fiber_diam']).agg(np.max)
print("Percent Error by Fiber Diameter")
print("Mean: ", pe_means['pe'])
print("Median: ", pe_medians['pe'])
print("Min: ", pe_mins['pe'])
print("Max: ", pe_maxs['pe'])
# now do the same but for absolute error new column 'ae' is the absolute value of column 'pe'
deffinalmatch['ae'] = deffinalmatch['pe'].abs()
ae_means = deffinalmatch.groupby(['fiber_diam']).agg(np.mean)
ae_medians = deffinalmatch.groupby(['fiber_diam']).agg(np.median)
ae_mins = deffinalmatch.groupby(['fiber_diam']).agg(np.min)
ae_maxs = deffinalmatch.groupby(['fiber_diam']).agg(np.max)
print("Absolute Error by Fiber Diameter")
print("Mean: ", ae_means['ae'])
print("Median: ", ae_medians['pe'])
print("Min: ", ae_mins['ae'])
print("Max: ", ae_maxs['ae'])
# %% r2 vals
thisdef = deftomatch.query('deformation!="ASCENT" and fiber_diam in [3,13]')
thisdef['deformed'] = thisdef['deformation'] == "Structural"

comparecols = [
    'peak_second_diff',
    'peri_thk',
    'peri_thk_act_site',
    'tortuosity',
    'minimum_efib_distance',
    'apnode_efib_distance',
]
# comparecols = ['peak_second_diff','peri_thk']

alldf = []

for name in comparecols:
    # Calculate correlation for each column with the "threshold" column based on the grouping variables
    correlation_df = (
        thisdef.groupby(["type", "deformed", "fiber_diam", "nerve_label"])[[name, "threshold"]].corr().reset_index()
    )
    correlation_df['name'] = name
    alldf.append(correlation_df.query('threshold!=1'))
alldf = pd.concat(alldf)
alldf.threshold = alldf.threshold**2
# Plot each column as a separate line in the line plot
g = sns.FacetGrid(data=alldf, col='type', row='fiber_diam', margin_titles=True)
g.map_dataframe(sns.pointplot, dodge=0.3, y='threshold', x='deformed', hue='name', palette='Set2')
leg = plt.legend(bbox_to_anchor=(1.4, -0.4))
new_labs = [
    'Peak Second Difference',
    'Perineurium Thickness (Center)',
    'Perineurium Thickness (Activation Site)',
    'Tortuosity',
    'Minimum Electrode-Fiber Distance',
    'AP Init Node Electrode-Fiber Distance',
]
for t, l in zip(leg.get_texts(), new_labs):
    t.set_text(l)
plt.ylim(0, 1)
g.set_ylabels(r'$R^2$', rotation=90)
g.set_titles(col_template="{col_name}", row_template="D: {row_name} μm")
for ax in g.axes[1, :]:
    ax.set_xticklabels(['No', 'Yes'])
    ax.set_xlabel("Deformed?")
# %% GLM
thisdef = deftomatch.query('deformation!="ASCENT" and fiber_diam in [3,13]')
thisdef['deformed'] = thisdef['deformation'] == "Structural"

comparecols = [
    'peak_second_diff',
    'peri_thk',
    'peri_thk_act_site',
    'tortuosity',
    'minimum_efib_distance',
    'apnode_efib_distance',
]
# comparecols = ['peak_second_diff','peri_thk']

alldf = []

for name in comparecols:
    # Calculate correlation for each column with the "threshold" column based on the grouping variables
    correlation_df = (
        thisdef.groupby(["type", "deformed", "fiber_diam", "nerve_label"])[[name, "threshold"]].corr().reset_index()
    )
    correlation_df['name'] = name
    alldf.append(correlation_df.query('threshold!=1'))
alldf = pd.concat(alldf)
alldf.threshold = alldf.threshold**2
# Plot each column as a separate line in the line plot
g = sns.FacetGrid(data=alldf, col='type', row='fiber_diam', margin_titles=True)
g.map_dataframe(sns.pointplot, dodge=0.3, y='threshold', x='deformed', hue='name', palette='Set2')
leg = plt.legend(bbox_to_anchor=(1.4, -0.4))
new_labs = [
    'Peak Second Difference',
    'Perineurium Thickness (Center)',
    'Perineurium Thickness (Activation Site)',
    'Tortuosity',
    'Minimum Electrode-Fiber Distance',
    'AP Init Node Electrode-Fiber Distance',
]
for t, l in zip(leg.get_texts(), new_labs):
    t.set_text(l)
plt.ylim(0, 1)
g.set_ylabels(r'$R^2$', rotation=90)
g.set_titles(col_template="{col_name}", row_template="D: {row_name} μm")
# %% r2 plots
sns.reset_orig()
thisdef = deftomatch.query('deformation!="ASCENT" and fiber_diam in [3,13]')
thisdef['deformed'] = thisdef['deformation'] == "Structural"

comparecols = [
    'peak_second_diff',
    'peak_abs_second_diff',
    'peri_thk',
    'peri_thk_act_site',
    'tortuosity',
    'minimum_efib_distance',
    'apnode_efib_distance',
]
# comparecols = ['peak_second_diff']

alldf = []

for name in comparecols:
    # plot threshold vs the comparecol
    g = sns.FacetGrid(data=thisdef, col='type', row='fiber_diam', margin_titles=True, sharey=False, sharex=False)
    g.map_dataframe(sns.scatterplot, y='threshold', x=name, hue='nerve_label', palette='Set2')
    # now plot with log scale x
    g = sns.FacetGrid(data=thisdef, col='type', row='fiber_diam', margin_titles=True, sharey=False, sharex=False)
    g.map_dataframe(sns.scatterplot, y='threshold', x=name, hue='deformed', palette='Set2')
    for ax in g.axes.flat:
        ax.set_xscale('log')
        ax.set_yscale('log')
    # do same R2 analysis as above, but take the log of both columns
    thisdef['log_threshold'] = np.log10(thisdef['threshold'])
    thisdef['log_' + name] = np.log10(thisdef[name])
    correlation_df = (
        thisdef.groupby(["type", "deformed", "fiber_diam", "nerve_label"])[['log_' + name, "log_threshold"]]
        .corr()
        .reset_index()
    )
    correlation_df['name'] = name
    alldf.append(correlation_df.query('log_threshold!=1'))
alldf = pd.concat(alldf)
alldf.log_threshold = alldf.log_threshold**2
# Plot each column as a separate line in the line plot
sns.set(font_scale=1.75, style='whitegrid')
g = sns.FacetGrid(data=alldf, col='type', row='fiber_diam', margin_titles=True)
g.map_dataframe(sns.pointplot, dodge=0.3, y='log_threshold', x='deformed', hue='name', palette='Set2')
leg = plt.legend(bbox_to_anchor=(1.4, -0.4))
new_labs = [
    'Peak Second Difference',
    'Peak Absolute Second Difference',
    'Perineurium Thickness (Center)',
    'Perineurium Thickness (Activation Site)',
    'Tortuosity',
    'Minimum Electrode-Fiber Distance',
    'AP Init Node Electrode-Fiber Distance',
]
for t, l in zip(leg.get_texts(), new_labs):
    t.set_text(l)
plt.ylim(0, 1)
g.set_ylabels(r'$R^2$', rotation=90)
g.set_titles(col_template="{col_name}", row_template="D: {row_name} μm")

# %% new idea analysis pick lowthresh
grouped = concats.query('contact != "center"').groupby(
    ['nerve_label', 'fiber_diam', 'master_fiber_index', 'deformation']
)
# take min threshold and threshold3d for each group
minthresh = grouped['threshold'].min().reset_index()
minthresh3d = grouped['threshold3d'].min().reset_index()
# merge minthresh and minthresh3d
minthresh.dropna(inplace=True)
minthresh = minthresh.merge(
    minthresh3d, on=['nerve_label', 'fiber_diam', 'master_fiber_index', 'deformation'], suffixes=['_2d', '_3d']
)
# plot minthresh vs minthresh3d
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

nsimdata = minthresh.query('fiber_diam in [3]')  # TODO replace all cath comparison with non
g = sns.relplot(
    data=nsimdata.rename(columns={'nerve_label': 'Sample'}),
    kind='scatter',
    col='deformation',
    x='threshold',
    y='threshold3d',
    # hue='Sample',
    color='white',
    # s=20,
    palette='colorblind',
    facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
    edgecolor='black',
    linewidth=1,
    alpha=1,
)
# TODO Clean up this calc
for diam, pos, ax, deformation in zip([3, 3], (0.2, 0.8, 0.2, 0.8), g.axes.ravel(), ["Undeformed", "Structural"]):
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    lim = np.amax([ax.get_xlim()[1], ax.get_ylim()[1]])
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} μm: {r ** 2:.2f}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='unity line')
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])

    # ax.set_yticks(ax.get_xticks())
g.set_titles('D: {col_name} μm')
g.set_xlabels('extrusion Threshold (mA)')
g.set_ylabels('true-3D Threshold (mA)')
g.set_titles(row_template='', col_template='Deformation: {col_name}')
# %% CCC comparison maincomp
nsimdata = concats.query(
    f'fiber_diam in [3,13] and contact in {main_comparison}'
)  # TODO replace all cath comparison with non
g = sns.relplot(
    data=nsimdata.rename(columns={'nerve_label': 'Sample'}),
    kind='scatter',
    col='fiber_diam',
    x='threshold',
    y='threshold3d',
    # hue='Sample',
    color='white',
    row='deformation',
    # s=20,
    palette='colorblind',
    facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
    edgecolor='black',
    linewidth=1,
    alpha=1,
)
# TODO Clean up this calc
for diam, pos, ax, deformation in zip(
    [3, 13, 3, 13], (0.2, 0.8, 0.2, 0.8), g.axes.ravel(), ["Undeformed", "Undeformed", "Structural", "Structural"]
):
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    lim = ax.get_xlim()[1]
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} μm: {r ** 2:.2f}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='unity line')
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])

    # ax.set_yticks(ax.get_xticks())
g.set_titles(col_template='D: {col_name} μm')
g.set_xlabels('extrusion Threshold (mA)')
g.set_ylabels('true-3D Threshold (mA)')
# %% remake the above with individual calls of histplot
sns.set(font_scale=1.75, style='whitegrid')
newthreshz = newdefdat.query("deformation=='Structural'")
newthreshz['activation_zpos_oneten'] = newthreshz['activation_zpos_oneten'] / 10000
fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
for nerve_label in pd.unique(newthreshz.nerve_label):
    for ax, modeltype in zip(axs, ["extrusion", "true-3D"]):
        g = sns.histplot(
            data=newthreshz.query(
                f"fiber_diam in [3,13] and contact in {main_comparison} and nerve_label=='{nerve_label}' and type=='{modeltype}'"
            ).rename(columns={'nerve_label': 'Sample'}),
            y='activation_zpos_oneten',
            hue='fiber_diam',
            # hue='Sample',
            # facet_kws={'sharex': False},
            # kind='kde',
            palette=[sns.color_palette('binary')[5], sns.color_palette('binary')[2]],
            common_norm=False,
            # legend=False,
            # multiple="fill",
            element='poly',
            fill=False,
            bins=np.arange(1.5, 3.6, 0.1),
            ax=ax,
        )
# delete both legends and remake my own
for ax in axs:
    ax.get_legend().remove()
# make my own legend
axs[0].plot([], [], color=sns.color_palette('binary')[5], label='3 μm', linewidth=2)
axs[0].plot([], [], color=sns.color_palette('binary')[2], label='13 μm', linewidth=2)
# put legenbd to right of figure
axs[0].axhspan(2.005, 2.205, color='red', alpha=0.2, label='contacts')
axs[1].axhspan(2.005, 2.205, color='red', alpha=0.2, label='_')
axs[0].axhspan(2.805, 3.005, color='red', alpha=0.2, label='_')
axs[1].axhspan(2.805, 3.005, color='red', alpha=0.2, label='_')

axs[0].legend(loc='center left', bbox_to_anchor=(2.2, 0.5))
axs[0].set_xlim(reversed(axs[0].get_xlim()))
axs[0].set_ylabel("Activation Location (cm)")
axs[0].set_title("extrusion-150%")
axs[1].set_title("true-3D-150%")
# %% strength duration correlation
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = threshload.query('contact=="3D"')
for comparison in [
    ['threshold', 'tortuosity'],
    ['threshold', 'nodal_tortuosity'],
    ['threshold', 'cuff_tortuosity'],
    ['threshold', 'peri_thk_act_site'],
    ['threshold', 'smallest_thk_under_cuff'],
    ['threshold', 'peri_thk_act_site'],
    ['smallest_thk_under_cuff', 'peri_thk_act_site'],
    # TODO try smallest thk over whole cuffspan as well as over both cuff spans
    ['threshold', 'minimum_efib_distance'],
    ['threshold', 'peak_second_diff'],
    ['peak_second_diff', 'peri_thk_act_site'],
    ['peak_second_diff', 'tortuosity'],
    ['peak_second_diff', 'nodal_tortuosity'],
    ['peak_second_z', 'activation_zpos'],
    ['peak_second_long_pos', 'long_ap_pos'],
    ['peak_second_diff_node', 'apnode'],
    ['peak_second_diff', 'minimum_efib_distance'],
]:
    corrs = (
        usedata.groupby(['sample', 'fiber_diam', 'pulse_width', 'contact', 'nerve_label'])[comparison]
        .corr()
        .iloc[0::2, -1]
    )
    corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
    corrs['pulse_width'] = pd.Categorical(corrs['pulse_width'].astype(float), ordered=True)
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    corrs['contact'] = pd.Categorical(corrs['contact'], ordered=True)
    plt.figure()
    # g = sns.FacetGrid(data=corrs,col='contact')
    # g.map_dataframe(sns.stripplot, hue='Sample', y='correlation',dodge=True, palette='colorblind')
    # g.map_dataframe(sns.boxplot,y='correlation', boxprops={'facecolor':'None'},whis=100)
    sns.swarmplot(data=corrs, y='fiber_diam', hue='pulse_width', x='correlation', dodge=True, s=6)
    sns.boxplot(data=corrs, y='fiber_diam', x='correlation', boxprops={'facecolor': 'None'}, whis=100)
    # plt.subplots_adjust(top=0.8)
    plt.title(f'Correlation between \n{comparison[0]} and {comparison[1]}', pad=25)
    plt.gca().set_xlim([-1, 1])
    plt.legend(title="PW (ms)", bbox_to_anchor=(1, 1))
    # plt.ylabel('')
    # plt.gca().set_yticklabels('')
    plt.gcf().set_size_inches([6, 5])
# %% strength duration dose-response
peses = []
pemeans = []
for stringdat in ['Undeformed', 'ASCENT', 'Structural']:
    subdat = newdefdat.query(f"deformation=='{stringdat}' and contact in {main_comparison}")
    sns.set(font_scale=1.75)
    sns.set_style('whitegrid')
    plt.figure()
    levels = {
        'onset': 10,
        'saturation': 90,
    }
    subdat['fiber_diam'] = subdat['fiber_diam'].astype(int)
    subdat['pulse_width'] = subdat['pulse_width'].astype(float)
    subdat['threshold'] = subdat['threshold'].astype(float)

    grouped = subdat.groupby(
        ['sample', 'fiber_diam', 'pulse_width', 'type', 'sim', 'nerve_label', 'model', 'nsim', 'deformation']
    )
    analysis = grouped.agg(
        {
            'threshold': [
                lambda x: np.percentile(x, q=levels['onset']),
                np.median,
                lambda x: np.percentile(x, q=levels['saturation']),
            ]
        }
    )
    analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
    analysis.rename(columns={'threshold_<lambda_0>': 'onset', 'threshold_<lambda_1>': 'saturation'}, inplace=True)
    analysis = analysis.reset_index()
    # combine onset, saturation, and half into one column with identifier
    compiled_data = analysis.melt(
        id_vars=['sample', 'fiber_diam', 'pulse_width', 'sim', 'type', 'nerve_label', 'model', 'nsim'],
        value_vars=['onset', 'saturation'],
        var_name='level',
        value_name='threshold',
    )

    # set up facetgrid with nsim as row and level as columns
    compiled_data.reset_index(inplace=True)
    # set fiber_diam to category
    compiled_data.type = compiled_data.type.astype('category')
    # remove all rows where modulus sample with 0 is 0
    compiled_data = compiled_data.query('sample % 10 != 0')
    # add a units column with unique number for each combination of fiber_diam and level
    compiled_data['units'] = compiled_data.groupby(['fiber_diam', 'level', 'nerve_label', 'pulse_width']).ngroup()
    compiled_data['fiber_diam'] = compiled_data['fiber_diam'].astype(int)
    compiled_data['pulse_width'] = compiled_data['pulse_width'].astype(float)

    compiled_data.dropna(inplace=True)

    # calculate percent error for sample onset and saturation as well as population onset and saturation
    pes = []
    for level in ['onset', 'saturation']:
        for sam in compiled_data.nerve_label.unique():
            for fiber_diam in compiled_data.fiber_diam.unique():
                for pulse_width in compiled_data.pulse_width.unique():
                    a = compiled_data.query(
                        f"nerve_label == '{sam}' and level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam} and pulse_width == {pulse_width}"
                    )['threshold'].values
                    b = compiled_data.query(
                        f"nerve_label == '{sam}' and level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam} and pulse_width == {pulse_width}"
                    )['threshold'].values
                    assert len(a) == len(b) == 1
                    pe_res = pe(a[0], b[0])
                    pes.append(
                        {
                            'level': level,
                            'nerve_label': sam,
                            'fiber_diam': fiber_diam,
                            'pulse_width': pulse_width,
                            'pe': pe_res,
                        }
                    )
    pes = pd.DataFrame(pes)
    pes["deformation"] = stringdat
    peses.append(pes)
    print('Max 3 um', stringdat, np.amax(pes.query('fiber_diam==3').pe))
    print('Max 13 um', stringdat, np.amax(pes.query('fiber_diam==13').pe))

    # # now calculate percent error for population onset and saturation
    # pemean = []
    # for level in ['onset', 'saturation']:
    #     for fiber_diam in compiled_data.fiber_diam.unique():
    #         a = compiled_data.query(f"level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam")[
    #             'threshold'
    #         ].values
    #         b = compiled_data.query(f"level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}")[
    #             'threshold'
    #         ].values
    #         pe_res = pe(np.median(a), np.median(b))
    #         pemean.append({'level': level, 'fiber_diam': fiber_diam, 'pe': pe_res})

    # pemean = pd.DataFrame(pemean)
    # pemean["deformation"] = stringdat
    # pemeans.append(pemean)
    # print('Max 3 um median', stringdat, np.amax(pemean.query('fiber_diam==3').pe))
    # print('Max 13 um median', stringdat, np.amax(pemean.query('fiber_diam==13').pe))

sns.set(font_scale=2, style='white')
allpes = pd.concat(peses)
g = sns.relplot(
    data=allpes,
    kind='line',
    col='deformation',
    style='level',
    y='pe',
    x='pulse_width',
    markers=True,
    color='black',
    units='nerve_label',
    row='fiber_diam',
    estimator=None,
    linewidth=2,
    markersize=8,
    # hue='nerve_label'
    facet_kws=dict(margin_titles=True),
)

g.set_xlabels('Pulse Width (ms)')
g.set_titles(col_template='Deformation: {col_name}')
g.set_ylabels("Percent Difference")
plt.ylim([0, None])
ylimpe = g.axes[0][0].get_ylim()

slotkin = (
    compiled_data.query('fiber_diam==3')
    .sort_values(['type', 'nerve_label'])
    .query('level=="onset"')
    .drop(columns=['sample', 'index', 'nsim', 'model', 'units'])
)
print(slotkin)
slotkin2 = (
    compiled_data.query('fiber_diam==3')
    .sort_values(['type', 'nerve_label'])
    .query('level=="saturation"')
    .drop(columns=['sample', 'index', 'nsim', 'model', 'units'])
)
print(slotkin)
# %% dose-response example
plt.figure()
drthis = defdr.query(
    f"fiber_diam in [3] and contact in {cath_comparison} and nerve_label =='6R' and deformation=='Structural'"
)

# Query the rows with type '3D'
df_3d = drthis.query("modeltype == 'true-3D'")

# Query the rows with type '2D' and sample ending in '1'
df_2d = drthis.query("modeltype == 'extrusion'")

# Merge the 3D and 2D data, keeping track of original row indices
merged_df = pd.merge(
    df_3d, df_2d, on=['nerve_label', 'master_fiber_index', 'nsim', 'deformation'], suffixes=('_3d', '_2d'), how='left'
)  # TODO remove this how

# Update the 'inner', 'outer', and 'fiber' columns in the original DataFrame
drthis.loc[df_3d.index, 'inner'] = merged_df['inner_2d'].values
drthis.loc[df_3d.index, 'outer'] = merged_df['outer_2d'].values
drthis.loc[df_3d.index, 'fiber'] = merged_df['fiber_2d'].values

sns.set(font_scale=1.25, style='white')
g = sns.scatterplot(
    data=drthis,
    y='percent_activated',
    x='threshold',
    hue='inner',
    palette='rainbow',
    linewidth=0,
    s=200,
    marker=r'$\backslash$',
    # legend=False
)

sns.lineplot(
    data=drthis,
    y='percent_activated',
    x='threshold',
    color='k',
    linewidth=2,
    estimator=None,
    style='modeltype',
    legend=False,
    zorder=1,
    alpha=1,
)


plt.xlabel('Threshold (mA)')
plt.ylim([0, 1])
# plt.gcf().set_size_inches(8,4)
plt.ylabel('Proportion Activated')
# create legend, circle = extrusion, X= true-3D
# create handles manually
from matplotlib.lines import Line2D

legend_elements = [
    Line2D([0], [0], label='extrusion', color='k', linestyle='-', alpha=0.6),
    Line2D([0], [0], label='true-3D', color='k', linestyle='--', alpha=0.6),
]
legend_labels = ['extrusion', 'true-3D']
g.legend(handles=legend_elements, labels=legend_labels, loc='lower right')
plt.figure()
g = sns.swarmplot(
    data=drthis,
    y='threshold',
    x='modeltype',
    hue='inner',
    palette='rainbow',
    linewidth=0,
    s=5,
    # legend=False
)
plt.xlabel('')
plt.ylabel('Threshold (mA)')
plt.legend([], [], frameon=False)
plt.ylim(0, None)
# %%
nowcomp = concats.query(
    f"fiber_diam in [3] and contact in {cath_comparison} and nerve_label =='6R' and deformation=='Structural'"
)
sns.scatterplot(nowcomp, x='threshold', y='threshold3d', hue='inner', palette='rainbow')
ax = plt.gca()
lim = np.amax([ax.get_xlim()[1], ax.get_ylim()[1]])
# add correlation to plot
# ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
# ax.set_title(f'{diam} μm')
ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='unity line')
plt.legend()
# ax.set_aspect('equal', 'box')
# ax.apply_aspect()
ax.set_xlim([0, lim])
ax.set_ylim([0, lim])
ax.set_aspect('equal')
plt.legend([], [], frameon=False)
plt.gcf().set_size_inches(5, 5)
plt.ylabel('True-3D Threshold (mA)')
plt.xlabel('Extrusion threshold (mA)')

# %% activation order fiberdiam
newdefdr = defdr.copy()
newdefdr = newdefdr.query("'6R' in nerve_label and deformation!='ASCENT'").sort_values('modeltype')
newdefdr['deformation'] = pd.Categorical(newdefdr['deformation'], categories=['Undeformed', 'Structural'])
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
newdefdr['percent_activated3'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newdefdr.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        f'fiber_diam == 3 and nerve_label == @row.nerve_label and modeltype == @row.modeltype and sim == @row.sim and master_fiber_index == @row.master_fiber_index and contact in {cath_comparison} and deformation==@row.deformation'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newdefdr.loc[row.Index, 'percent_activated3'] = val
# %% plot
sns.set(font_scale=1.5, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    # row='deformation',
    data=newdefdr.query(f"fiber_diam in [3,13] and contact in {cath_comparison} and deformation=='Structural'"),
    y='percent_activated',
    x='fiber_diam',
    units='nerve_label',
    col='modeltype',
    palette='plasma',
    hue='percent_activated3',
    # hue='inner',
    estimator=None,
    linewidth=0,
    facet_kws={'margin_titles': True},
    s=25,
)
plt.subplots_adjust(top=0.87)
# plt.suptitle(stringdat, x=0.37)
g.set_titles(row_template='', col_template='{col_name}')
g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion fibers activated')
g.set_xlabels('D (μm)')
g.legend.remove()
norm = plt.Normalize(0, 1)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])

# Remove the legend and add a colorbar
g.figure.colorbar(
    sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.8, label='Proportion 3 μm\nfibers activated', pad=0.1
).ax.yaxis.set_ticks_position('left')

for i, mod in enumerate(newdefdr.modeltype.sort_values().unique()):
    shortdat = newdefdr.query(
        f"modeltype==@mod and fiber_diam in [13] and contact in {cath_comparison} and deformation=='Structural'"
    )
    data2d = shortdat.sort_values('percent_activated').master_fiber_index
    data3d = shortdat.sort_values('percent_activated3').master_fiber_index
    rc = compute_reorder_cost(list(data2d), list(data3d))
    ax = g.axes[0][i]
    text = ax.get_title() + f' - AR: {round(rc,3)}'
    ax.set_title(text)

# %% CCC comparison maincomp
nsimdata = concats.query(
    f'fiber_diam in [3] and contact in {main_comparison}'
)  # TODO replace all cath comparison with non
g = sns.relplot(
    data=nsimdata.rename(columns={'nerve_label': 'Sample'}),
    kind='scatter',
    col='deformation',
    x='threshold',
    y='threshold3d',
    # hue='Sample',
    color='white',
    # s=20,
    palette='colorblind',
    facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
    edgecolor='black',
    linewidth=1,
    alpha=1,
)
lim = {}
# TODO Clean up this calc
for diam, pos, ax, deformation in zip([3, 3], (0.2, 0.8, 0.2, 0.8), g.axes.ravel(), ["Undeformed", "Structural"]):
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    lim[deformation] = np.amax([ax.get_xlim()[1], ax.get_ylim()[1]])
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} μm: {r ** 2:.2f}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim[deformation]], [0, lim[deformation]], '--k', linewidth=2, label='unity line')
    plt.legend()
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim[deformation]])
    ax.set_ylim([0, lim[deformation]])

    # ax.set_yticks(ax.get_xticks())
g.set_titles('D: {col_name} μm')
g.set_xlabels('extrusion Threshold (mA)')
g.set_ylabels('true-3D Threshold (mA)')
g.set_titles(row_template='', col_template='Deformation: {col_name}')

# new idea analysis pick lowthresh
grouped = concats.query('contact != "center"').groupby(
    ['nerve_label', 'fiber_diam', 'master_fiber_index', 'deformation']
)
# take min threshold and threshold3d for each group
minthresh = grouped['threshold'].min().reset_index()
minthresh3d = grouped['threshold3d'].min().reset_index()
# merge minthresh and minthresh3d
minthresh.dropna(inplace=True)
minthresh = minthresh.merge(
    minthresh3d, on=['nerve_label', 'fiber_diam', 'master_fiber_index', 'deformation'], suffixes=['_2d', '_3d']
)
# plot minthresh vs minthresh3d
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

nsimdata = minthresh.query('fiber_diam in [3]')  # TODO replace all cath comparison with non
g = sns.relplot(
    data=nsimdata.rename(columns={'nerve_label': 'Sample'}),
    kind='scatter',
    col='deformation',
    x='threshold',
    y='threshold3d',
    # hue='Sample',
    color='white',
    # s=20,
    palette='colorblind',
    facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
    edgecolor='black',
    linewidth=1,
    alpha=1,
)
# TODO Clean up this calc
for diam, pos, ax, deformation in zip([3, 3], (0.2, 0.8, 0.2, 0.8), g.axes.ravel(), ["Undeformed", "Structural"]):
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} μm: {r ** 2:.2f}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim[deformation]], [0, lim[deformation]], '--k', linewidth=2, label='unity line')
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim[deformation]])
    ax.set_ylim([0, lim[deformation]])
    plt.legend()
    # ax.set_yticks(ax.get_xticks())
g.set_titles('D: {col_name} μm')
g.set_xlabels('extrusion Threshold (mA)')
g.set_ylabels('true-3D Threshold (mA)')
g.set_titles(row_template='', col_template='Deformation: {col_name}')
# %% dose-response slotkin
peses = []
pemeans = []
subdat = threshload.query(f"contact in {main_comparison}")
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
plt.figure()
levels = {
    'onset': 10,
    'saturation': 90,
}
subdat['fiber_diam'] = subdat['fiber_diam'].astype(int)
subdat['pulse_width'] = subdat['pulse_width'].astype(float)
subdat['threshold'] = subdat['threshold'].astype(float)

grouped = subdat.groupby(['sample', 'fiber_diam', 'pulse_width', 'type', 'sim', 'nerve_label', 'model', 'nsim'])
analysis = grouped.agg(
    {
        'threshold': [
            lambda x: np.percentile(x, q=levels['onset']),
            np.median,
            lambda x: np.percentile(x, q=levels['saturation']),
        ]
    }
)
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.rename(columns={'threshold_<lambda_0>': 'onset', 'threshold_<lambda_1>': 'saturation'}, inplace=True)
analysis = analysis.reset_index()
# combine onset, saturation, and half into one column with identifier
compiled_data = analysis.melt(
    id_vars=['sample', 'fiber_diam', 'pulse_width', 'sim', 'type', 'nerve_label', 'model', 'nsim'],
    value_vars=['onset', 'saturation'],
    var_name='level',
    value_name='threshold',
)

# set up facetgrid with nsim as row and level as columns
compiled_data.reset_index(inplace=True)
# set fiber_diam to category
compiled_data.type = compiled_data.type.astype('category')
# remove all rows where modulus sample with 0 is 0
compiled_data = compiled_data.query('sample % 10 != 0')
# add a units column with unique number for each combination of fiber_diam and level
compiled_data['units'] = compiled_data.groupby(['fiber_diam', 'level', 'nerve_label', 'pulse_width']).ngroup()
compiled_data['fiber_diam'] = compiled_data['fiber_diam'].astype(int)
compiled_data['pulse_width'] = compiled_data['pulse_width'].astype(float)

compiled_data.dropna(inplace=True)

# calculate percent error for sample onset and saturation as well as population onset and saturation
pes = []
for level in ['onset', 'saturation']:
    for sam in compiled_data.nerve_label.unique():
        for fiber_diam in compiled_data.fiber_diam.unique():
            for pulse_width in compiled_data.pulse_width.unique():
                a = compiled_data.query(
                    f"nerve_label == '{sam}' and level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam} and pulse_width == {pulse_width}"
                )['threshold'].values
                b = compiled_data.query(
                    f"nerve_label == '{sam}' and level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam} and pulse_width == {pulse_width}"
                )['threshold'].values
                assert len(a) == len(b) == 1
                pe_res = pe(a[0], b[0])
                pes.append(
                    {
                        'level': level,
                        'nerve_label': sam,
                        'fiber_diam': fiber_diam,
                        'pulse_width': pulse_width,
                        'pe': pe_res,
                    }
                )
pes = pd.DataFrame(pes)
pes["deformation"] = stringdat
peses.append(pes)
print('Max 3 um', stringdat, np.amax(pes.query('fiber_diam==3').pe))
print('Max 13 um', stringdat, np.amax(pes.query('fiber_diam==13').pe))

# # now calculate percent error for population onset and saturation
# pemean = []
# for level in ['onset', 'saturation']:
#     for fiber_diam in compiled_data.fiber_diam.unique():
#         a = compiled_data.query(f"level == '{level}' and type == 'extrusion' and fiber_diam == {fiber_diam")[
#             'threshold'
#         ].values
#         b = compiled_data.query(f"level == '{level}' and type == 'true-3D' and fiber_diam == {fiber_diam}")[
#             'threshold'
#         ].values
#         pe_res = pe(np.median(a), np.median(b))
#         pemean.append({'level': level, 'fiber_diam': fiber_diam, 'pe': pe_res})

# pemean = pd.DataFrame(pemean)
# pemean["deformation"] = stringdat
# pemeans.append(pemean)
# print('Max 3 um median', stringdat, np.amax(pemean.query('fiber_diam==3').pe))
# print('Max 13 um median', stringdat, np.amax(pemean.query('fiber_diam==13').pe))

sns.set(font_scale=2, style='white')
allpes = pd.concat(peses)
g = sns.relplot(
    data=allpes,
    kind='line',
    col='deformation',
    style='level',
    y='pe',
    x='pulse_width',
    markers=True,
    color='black',
    units='nerve_label',
    row='fiber_diam',
    estimator=None,
    linewidth=2,
    markersize=8,
    # hue='nerve_label'
    facet_kws=dict(margin_titles=True),
)

g.set_xlabels('Pulse Width (ms)')
g.set_titles(col_template='Deformation: {col_name}')
g.set_ylabels("Percent Difference")
plt.ylim([0, None])
ylimpe = g.axes[0][0].get_ylim()

slotkin = (
    compiled_data.query('fiber_diam==3')
    .sort_values(['type', 'nerve_label'])
    .query('level=="onset"')
    .drop(columns=['sample', 'index', 'nsim', 'model', 'units'])
)
print(slotkin)
slotkin2 = (
    compiled_data.query('fiber_diam==3')
    .sort_values(['type', 'nerve_label'])
    .query('level=="saturation"')
    .drop(columns=['sample', 'index', 'nsim', 'model', 'units'])
)
print(slotkin)
# %% activation residual
# TODO fix
import seaborn as sns

sns.set(font_scale=1.5, style='whitegrid')
tcopy = threshload.copy()
tcopy['type'] = pd.Categorical(tcopy['type'], categories=['true-3D', 'extrusion'], ordered=True)
tcopy['zres'] = (threshload.activation_zpos - (50100 - threshload.peak_second_z)) / (threshload.fiber_diam * 100)
# g= sns.stripplot(data=tcopy,y='zres',x='type', hue='fiber_diam',dodge=True,palette='RdPu',linewidth=0,edgecolor='k')
# g= sns.lineplot(data=tcopy,y='zres',hue='type', x='fiber_diam',palette=pal2d3d,errorbar=('pi',100),marker='o')
g = sns.pointplot(
    data=tcopy,
    y='zres',
    hue='type',
    x='fiber_diam',
    palette=reversed(pal2d3d),
    errorbar=('pi', 100),
    dodge=0.1,
    estimator='median',
)
legend = plt.legend(title='', ncols=2)
# plt.xticks(threshload.fiber_diam.unique())
plt.xlabel('D (μm)')
plt.ylabel('Activation residual\n(no. node distances)')
# %% activation residual nodes
# TODO Fix
import seaborn as sns

sns.set(font_scale=1.5, style='whitegrid')
tcopy = threshload.copy()
tcopy['type'] = pd.Categorical(tcopy['type'], categories=['true-3D', 'extrusion'], ordered=True)
tcopy['zres'] = threshload.apnode - threshload.peak_second_diff_node
# g= sns.stripplot(data=tcopy,y='zres',x='type', hue='fiber_diam',dodge=True,palette='RdPu',linewidth=0,edgecolor='k')
# g= sns.lineplot(data=tcopy,y='zres',hue='type', x='fiber_diam',palette=pal2d3d,errorbar=('pi',100),marker='o')
g = sns.pointplot(
    data=tcopy,
    y='zres',
    hue='type',
    x='fiber_diam',
    palette=reversed(pal2d3d),
    errorbar=('pi', 100),
    dodge=0.1,
    estimator='median',
)
legend = plt.legend(title='', ncols=2)
# plt.xticks(threshload.fiber_diam.unique())
plt.xlabel('D (μm)')
plt.ylabel('Activation residual\n(no. node distances)')
# %% GLM

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

# Create GENERALIZED LINEAR MODEL
# formula: Default stands for intercept (dependet vriable)
# The intercept is the predicted value of the dependent variables, when all the independent variables /Risk level, YOB, DJX, GDP/ are 0.
# You can try family Tweedie with its link function log as well in this case
# ['peak_second_diff', 'peri_thk', 'peri_thk_act_site', 'tortuosity', 'minimum_efib_distance','apnode_efib_distance']

for deformation in ["Structural", "Undeformed"]:
    for modeltype in ['extrusion', 'true-3D']:
        for sample in pd.unique(thisdef['nerve_label']):

            print(deformation, modeltype, sample)
            model = smf.glm(
                formula="threshold ~ " + ' + '.join(comparecols),
                data=thisdef.query('deformation==@deformation and type==@modeltype and nerve_label==@sample'),
                family=sm.families.Gaussian(),
            )

            # Fit the model
            result = model.fit()
            # Display and interpret results
            print(result.summary())
            # Estimated default probabilities
            predictions = result.predict()
# %% organizatino compare type
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshdat = newdefdat.query(f'contact in {cath_comparison} and deformation=="Structural"')

scores = []
for nerve in pd.unique(threshdat['nerve_label']):
    for model in ['extrusion', 'true-3D']:
        shortdat = threshdat.query(f'nerve_label=="{nerve}" and type=="{model}"')
        datasmol = shortdat.query('nsim==0').sort_values('threshold').master_fiber_index
        databeeg = shortdat.query('nsim==5').sort_values('threshold').master_fiber_index
        rc = compute_reorder_cost(list(datasmol), list(databeeg))
        scores.append({'sample': nerve, 'type': model, 'scoresmolbeeg': rc})
fig = plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.stripplot(
    data=scoredat, y='scoresmolbeeg', x='type', color='black', size=10, jitter=False, hue='sample', palette=defpal
)
ax = sns.lineplot(data=scoredat, y='scoresmolbeeg', x='type', hue='sample', palette=defpal)
plt.ylabel('Activation Reordering')
plt.xlabel('')
plt.ylim(0, 0.6)
print(scoredat.groupby(['type']).median())
plt.xlim(-1, 2)
# plt.gcf().set_size_inches([3, 5])
fig.set_size_inches(4, 4)
plt.legend([], [], frameon=False)
plt.xticks(rotation=45)
# %% organizatino compare type and def
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshdat = newdefdat.query(f'contact in {cath_comparison}')

scores = []
for nerve in pd.unique(threshdat['nerve_label']):
    for model in ['extrusion', 'true-3D']:
        for deformtype in ['Undeformed', 'Structural']:
            shortdat = threshdat.query(f'nerve_label=="{nerve}" and type=="{model}" and deformation ==@deformtype')
            datasmol = shortdat.query('nsim==0').sort_values('threshold').master_fiber_index
            databeeg = shortdat.query('nsim==5').sort_values('threshold').master_fiber_index
            rc = compute_reorder_cost(list(datasmol), list(databeeg))
            scores.append({'sample': nerve, 'type': model, 'scoresmolbeeg': rc, 'deformation': deformtype})
fig = plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.stripplot(
    data=scoredat, y='scoresmolbeeg', x='type', color='black', size=10, jitter=False, hue='sample', palette=defpal
)
ax = sns.lineplot(data=scoredat, y='scoresmolbeeg', x='type', hue='sample', palette=defpal, style='deformation')
plt.ylabel('Activation Reordering')
plt.xlabel('')
plt.ylim(0, 0.6)
print(scoredat.groupby(['type']).median())
plt.xlim(-1, 2)
# plt.gcf().set_size_inches([3, 5])
fig.set_size_inches(4, 4)
plt.legend([], [], frameon=False)
plt.xticks(rotation=45)
ARdat = pd.DataFrame(scores)
# %% organizatino compare type and def
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshdat = newdefdat.query(f'contact in {cath_comparison}')

scores = []
for nerve in pd.unique(threshdat['nerve_label']):
    for model in ['extrusion', 'true-3D']:
        for deformtype in ['Undeformed', 'Structural']:
            shortdat = threshdat.query(f'nerve_label=="{nerve}" and type=="{model}" and deformation ==@deformtype')
            datasmol = shortdat.query('nsim==0').sort_values('threshold').master_fiber_index
            databeeg = shortdat.query('nsim==5').sort_values('threshold').master_fiber_index
            rc = compute_reorder_cost(list(datasmol), list(databeeg))
            scores.append({'sample': nerve, 'type': model, 'scoresmolbeeg': rc, 'deformation': deformtype})
fig = plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.stripplot(
    data=scoredat, y='scoresmolbeeg', x='type', color='black', size=10, jitter=False, hue='sample', palette=defpal
)
ax = sns.lineplot(data=scoredat, y='scoresmolbeeg', x='type', hue='sample', palette=defpal, style='deformation')
plt.ylabel('Activation Reordering')
plt.xlabel('')
plt.ylim(0, 0.6)
print(scoredat.groupby(['type']).median())
plt.xlim(-1, 2)
# plt.gcf().set_size_inches([3, 5])
fig.set_size_inches(4, 4)
plt.legend([], [], frameon=False)
plt.xticks(rotation=45)
ARdat = pd.DataFrame(scores)
# %% organizatino compare type and def
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshdat = newdefdat.copy()

scores = []
for nerve in pd.unique(threshdat['nerve_label']):
    for deformtype in ['Undeformed', 'Structural']:
        for contact in ['3D', 'cathodic', 'center', 'anodic']:
            shortdat = threshdat.query(f'nerve_label=="{nerve}" and deformation ==@deformtype and contact in @contact')
            datasmol = shortdat.query('nsim==0').sort_values('threshold').master_fiber_index
            databeeg = shortdat.query('nsim==5').sort_values('threshold').master_fiber_index
            rc = compute_reorder_cost(list(datasmol), list(databeeg))
            scores.append(
                {'sample': nerve, 'type': model, 'scoresmolbeeg': rc, 'deformation': deformtype, 'slice': contact}
            )
fig = plt.figure()
scoredat = pd.DataFrame(scores)
sns.stripplot(
    data=scoredat,
    y='scoresmolbeeg',
    x='slice',
    color='black',
    hue='deformation',
    dodge=True,
    linewidth=0.5,
    edgecolor='white',
)
sns.barplot(data=scoredat, y='scoresmolbeeg', x='slice', hue='deformation', palette='binary')
plt.ylabel('Activation Reordering')
plt.ylim(0, 0.6)
print(scoredat.groupby(['type']).median())
# plt.xlim(-1,2)
# plt.gcf().set_size_inches([3, 5])
plt.legend([], [], frameon=False)
plt.xticks(rotation=45)
plt.xlabel('extrusion slice')
g.set_titles(col_template='{col_name}')
g.fig.set_size_inches(8, 6)
# %% Percent Error deformed
sns.reset_orig()
sns.set(font_scale=1.5, style='white')
# apply pe to all rows of dataframe matched, with threshold3d as the correct value and threshold as the estimated value
deffinalmatch['pe'] = deffinalmatch.apply(lambda row: pe(row['threshold3d'], row['threshold']), axis=1)
plt.figure()
sns.barplot(data=deffinalmatch, x='nerve_label', y='pe', hue='fiber_diam', errorbar='se', palette="RdPu")
# plt.title('Threshold Percent Error by sample and fiber diameter')
legend = plt.legend(title='D (μm)', ncols=2, prop={'size': 12})
plt.xlabel('')
plt.ylabel('Percent Difference (%)')
plt.gcf().set_size_inches([6, 5])
legend.set_title('D (μm)', prop={'size': 12})

# calculate min, max, and mean percent error for each fiber diameter
pe_means = matched.groupby(['fiber_diam']).agg(np.mean)
pe_medians = matched.groupby(['fiber_diam']).agg(np.median)
pe_mins = matched.groupby(['fiber_diam']).agg(np.min)
pe_maxs = matched.groupby(['fiber_diam']).agg(np.max)
print("Percent Error by Fiber Diameter")
print("Mean: ", pe_means['pe'])
print("Median: ", pe_medians['pe'])
print("Min: ", pe_mins['pe'])
print("Max: ", pe_maxs['pe'])
# now do the same but for absolute error new column 'ae' is the absolute value of column 'pe'
matched['ae'] = matched['pe'].abs()
ae_means = matched.groupby(['fiber_diam']).agg(np.mean)
ae_medians = matched.groupby(['fiber_diam']).agg(np.median)
ae_mins = matched.groupby(['fiber_diam']).agg(np.min)
ae_maxs = matched.groupby(['fiber_diam']).agg(np.max)
print("Absolute Error by Fiber Diameter")
print("Mean: ", ae_means['ae'])
print("Median: ", ae_medians['pe'])
print("Min: ", ae_mins['ae'])
print("Max: ", ae_maxs['ae'])
plt.figure()
outdata = (
    deffinalmatch.query(f'contact in {main_comparison}')
    .groupby(['nerve_label', 'fiber_diam'])
    .agg({'threshold': ['mean', 'std', 'median'], 'threshold3d': ['mean', 'std', 'median']})
)
outdata.columns = ["_".join(col_name).rstrip('_') for col_name in outdata.columns]
outdata.reset_index(inplace=True)
outdata.dropna(inplace=True)
# calculate percent difference for mean, std, and median
outdata['mean_diff'] = (outdata.threshold_mean - outdata.threshold3d_mean) / outdata.threshold3d_mean * 100
outdata['std_diff'] = (outdata.threshold_std - outdata.threshold3d_std) / outdata.threshold3d_std * 100
outdata['median_diff'] = (outdata.threshold_median - outdata.threshold3d_median) / outdata.threshold3d_median * 100
# plot percent difference
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
sns.lineplot(data=outdata, x='fiber_diam', y='mean_diff', hue='nerve_label', palette='colorblind', marker='o')
plt.xticks(outdata.fiber_diam.unique())
plt.ylabel('Percent Difference (%)')
plt.xlabel('Fiber Diameter (μm)')
# plt.title('Mean Threshold Percent Difference')
plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0))
plt.show()
# %% CCC comparison maincomp
diam = 3
for comp in [an_comparison, main_comparison, cath_comparison]:
    deformation = "Structural"
    nsimdata = concats.query(
        f'fiber_diam in [3] and contact in {comp} and deformation==@deformation'
    )  # TODO replace all cath comparison with non
    g = sns.relplot(
        data=nsimdata.rename(columns={'nerve_label': 'Sample'}),
        kind='scatter',
        x='threshold',
        y='threshold3d',
        # hue='Sample',
        color='white',
        # s=20,
        facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
        edgecolor='black',
        linewidth=1,
        alpha=1,
    )
    lim = {}
    ax = g.axes[0, 0]
    # TODO Clean up this calc
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    lim[deformation] = np.amax([ax.get_xlim()[1], ax.get_ylim()[1]])
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} {comp[0]} μm CCC: {r ** 2:.2f}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim[deformation]], [0, lim[deformation]], '--k', linewidth=2, label='unity line')
    plt.legend()
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim[deformation]])
    ax.set_ylim([0, lim[deformation]])

    # ax.set_yticks(ax.get_xticks())
    g.set_titles('D: {col_name} μm')
    g.set_xlabels('extrusion Threshold (mA)')
    g.set_ylabels('true-3D Threshold (mA)')
    g.set_titles(row_template='', col_template='Deformation: {col_name}')

    g.set_titles(row_template='', col_template='Deformation: {col_name}')
    mid = [np.diff(plt.xlim()) / 2, np.diff(plt.ylim()) / 2]
    mid = [float(x) for x in mid]
    plt.arrow(mid[0] - 0.25, mid[1] + 0.25, -0.5, 0.5, color='black', width=0.08)
    plt.text(mid[0] - 1.5, mid[1] + 1, 'true-3D higher')
    plt.arrow(mid[0] + 0.25, mid[1] - 0.25, +0.5, -0.5, color='black', width=0.08)
    plt.text(mid[0], mid[1] - 1.25, 'extrusion higher')
# new idea analysis pick lowthresh
grouped = concats.query('contact != "center"').groupby(
    ['nerve_label', 'fiber_diam', 'master_fiber_index', 'deformation']
)
# take min threshold and threshold3d for each group
minthresh = grouped['threshold'].min().reset_index()
minthresh3d = grouped['threshold3d'].min().reset_index()
# merge minthresh and minthresh3d
minthresh.dropna(inplace=True)
minthresh = minthresh.merge(
    minthresh3d, on=['nerve_label', 'fiber_diam', 'master_fiber_index', 'deformation'], suffixes=['_2d', '_3d']
)
# plot minthresh vs minthresh3d
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

nsimdata = concats.query('fiber_diam in [3]')  # TODO replace all cath comparison with non
g = sns.relplot(
    data=nsimdata.rename(columns={'nerve_label': 'Sample'}).query('deformation==@deformation'),
    kind='scatter',
    x='threshold',
    y='threshold3d',
    # hue='Sample',
    color='white',
    # s=20,
    facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
    edgecolor='black',
    linewidth=1,
    alpha=1,
)
ax = g.axes[0, 0]
# TODO Clean up this calc
rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
# add correlation to plot
# ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
print(f'{diam} {deformation} μm CCC: {r ** 2:.2f}')
# ax.set_title(f'{diam} μm')
ax.plot([0, lim[deformation]], [0, lim[deformation]], '--k', linewidth=2, label='unity line')
# ax.set_aspect('equal', 'box')
# ax.apply_aspect()
ax.set_xlim([0, lim[deformation]])
ax.set_ylim([0, lim[deformation]])
plt.legend()
# ax.set_yticks(ax.get_xticks())
g.set_titles('D: {col_name} μm')
g.set_xlabels('extrusion Threshold (mA)')
g.set_ylabels('true-3D Threshold (mA)')
# %% activation order calc minthresh
newdefdr = defdr.query("nerve_label=='5R' and fiber_diam==3")
anode = newdefdr.query('contact=="anodic" and modeltype=="extrusion"')
cathode = newdefdr.query('contact=="cathodic" and modeltype=="extrusion"')
match_min = datamatch_merge(
    anode,
    cathode,
    'threshold',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index", "deformation"],
).drop(columns='contact')
match_min['threshold'] = np.minimum(match_min['threshold'], match_min['threshold3d'])
matchmindr = calculate_dose_response(
    match_min.drop(columns='percent_activated'),
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim'],
)
matchmindr['contact'] = 'minthresh'
newdefdr = pd.concat([newdefdr, matchmindr])
newdefdr = newdefdr.query("'5R' in nerve_label and deformation=='Structural'").sort_values('modeltype')
newdefdr['nerve_label'] = pd.Categorical(newdefdr['nerve_label'], categories=['5R'])
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
newdefdr['percent_activated2d'] = np.nan
newdefdr.reset_index(drop=True, inplace=True)
# go through every row and for each fiber find the 2D activation percent
for row in newdefdr.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        'modeltype == "true-3D" and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and master_fiber_index == @row.master_fiber_index and deformation==@row.deformation'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newdefdr.loc[row.Index, 'percent_activated3d'] = val
newdefdr['contact'] = pd.Categorical(
    newdefdr['contact'], categories=["3D", "cathodic", "center", "anodic", "minthresh"], ordered=True
)
# %% plot activation order minthresh
sns.set(font_scale=1.25, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    row='fiber_diam',
    data=newdefdr.query(f"fiber_diam in [3] and nerve_label=='5R'"),
    y='percent_activated',
    x='contact',
    units='nerve_label',
    palette='plasma',
    hue='percent_activated3d',
    # hue='inner',
    estimator=None,
    linewidth=0,
    facet_kws={'margin_titles': True},
    s=25,
)
plt.subplots_adjust(top=0.87)
# plt.suptitle(stringdat, x=0.37)
g.set_titles(row_template='', col_template='{col_name}')
g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion fibers activated')
g.set_xlabels('')
g.legend.remove()
norm = plt.Normalize(0, 1)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])
plt.axvline(0.5, linestyle='--', color='black', alpha=0.5)
# Remove the legend and add a colorbar
g.figure.colorbar(
    sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.8, label='Proportion true-3D\nfibers activated', pad=0.1
).ax.yaxis.set_ticks_position('left')
g.fig.set_size_inches(14, 4)
for i, con in enumerate(newdefdr.contact.sort_values().unique()[1:]):
    shortdat = newdefdr.query('fiber_diam==3 and contact==@con')
    data2d = shortdat.sort_values('percent_activated').master_fiber_index
    data3d = shortdat.sort_values('percent_activated3d').master_fiber_index
    rc = compute_reorder_cost(list(data2d), list(data3d))
    g.axes[0][0].text(i + 0.77, 1.065, round(rc, 3))
g.axes[0][0].text(-0.5, 1.07, 'Activation\nreordering: ')
# %% moreplots
sns.relplot(
    data=newdefdr, kind='scatter', x='percent_activated', y='percent_activated3d', hue='nerve_label', row='contact'
)

# %% calculate CCC for all comparisons, deformation types
CCC_vals = []
for deformation in ["Undeformed", "Structural"]:
    for diam in concats.fiber_diam.unique():
        for comp in [an_comparison, main_comparison, cath_comparison]:
            CCC_dat = concats.query('contact in @comp and deformation==@deformation and fiber_diam==@diam')
            r = concordance_correlation_coefficient(CCC_dat.threshold3d, CCC_dat.threshold)
            CCC_vals.append({'deformation': deformation, 'comparison': comp[0], 'CCC': r**2, 'diam': diam})
        # also do minthresh, where comp = "combo"
        CCC_dat = minthresh.query('deformation==@deformation and fiber_diam==@diam')
        r = concordance_correlation_coefficient(CCC_dat.threshold3d, CCC_dat.threshold)
        CCC_vals.append({'deformation': deformation, 'comparison': 'combo', 'CCC': r**2, 'diam': diam})
CCC_vals = pd.DataFrame(CCC_vals)
# plot
sns.set(font_scale=1.75)
sns.set_style('white')
g = sns.relplot(
    kind='line',
    data=CCC_vals,
    x='comparison',
    y='CCC',
    style='deformation',
    hue='diam',
    # palette='RdPu'
)
# %% newplot diff
sns.set(font_scale=1.5, style='ticks')
plt.figure()
outdata = (
    deffinalmatch.query(f'contact in {main_comparison}')
    .groupby(['nerve_label', 'fiber_diam'])
    .agg({'threshold': ['mean', 'std', 'median'], 'threshold3d': ['mean', 'std', 'median']})
)
outdata.columns = ["_".join(col_name).rstrip('_') for col_name in outdata.columns]
outdata.reset_index(inplace=True)
outdata.dropna(inplace=True)
# calculate percent difference for mean, std, and median
outdata['norm2d'] = outdata.threshold_mean / outdata.threshold3d_mean
# plot percent difference
sns.lineplot(data=outdata, x='fiber_diam', y='norm2d', hue='nerve_label', palette='colorblind', marker='o')
plt.xticks(outdata.fiber_diam.unique())
plt.ylabel('Normalized Mean\nextrusion Threshold')
plt.xlabel('Fiber Diameter (μm)')
plt.axhline(1, linestyle=':')
plt.ylim(0.6, 1.4)
# plt.title('Mean Threshold Percent Difference')
plt.legend(ncols=2)
plt.show()
# %% remake the above with individual calls of histplot
sns.set(font_scale=1.75, style='whitegrid')
newthreshz = newdefdat.query("deformation=='Structural'")
newthreshz['activation_zpos'] = 5.01 - newthreshz['peak_second_z'] / 10000
fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
for nerve_label in pd.unique(newthreshz.nerve_label):
    for ax, modeltype in zip(axs, ["extrusion", "true-3D"]):
        g = sns.histplot(
            data=newthreshz.query(
                f"fiber_diam in [3,13] and contact in {main_comparison} and nerve_label=='{nerve_label}' and type=='{modeltype}'"
            ).rename(columns={'nerve_label': 'Sample'}),
            y='activation_zpos',
            hue='fiber_diam',
            # hue='Sample',
            # facet_kws={'sharex': False},
            # kind='kde',
            palette=[sns.color_palette('binary')[5], sns.color_palette('binary')[2]],
            common_norm=False,
            # legend=False,
            # multiple="fill",
            element='poly',
            fill=False,
            bins=np.arange(1.5, 3.6, 0.1),
            ax=ax,
        )
# delete both legends and remake my own
for ax in axs:
    ax.get_legend().remove()
# make my own legend
axs[0].plot([], [], color=sns.color_palette('binary')[5], label='3 μm', linewidth=2)
axs[0].plot([], [], color=sns.color_palette('binary')[2], label='13 μm', linewidth=2)
# put legenbd to right of figure
axs[0].axhspan(2.005, 2.205, color='red', alpha=0.2, label='contacts')
axs[1].axhspan(2.005, 2.205, color='red', alpha=0.2, label='_')
axs[0].axhspan(2.805, 3.005, color='red', alpha=0.2, label='_')
axs[1].axhspan(2.805, 3.005, color='red', alpha=0.2, label='_')

axs[0].legend(loc='center left', bbox_to_anchor=(2.2, 0.5))
axs[0].set_xlim(reversed(axs[0].get_xlim()))
axs[0].set_ylabel("Peak Second 2nd diff. (cm)")
axs[0].set_title("extrusion")
axs[1].set_title("true-3D")
# %% Correlationall
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = threshload.copy()
for comparison in [
    # ['thresholdextrusion', 'threshold3d'],
    # ['thresholdextrusion', 'peri_thk'],
    # ['threshold3d', 'peri_thk'],
    # ['thresholdextrusion', 'minimum_efib_distance'],
    ['threshold', 'peak_second_diff'],
    # ['thresholdextrusion', 'peak_abs_second_diff'],
    # ['peak_second_diff', 'peri_thk'],
    # ['peak_second_diff', 'minimum_efib_distance'],
    # ['peak_second_z', 'activation_zpos'],
    # ['peak_second_diff_node', 'apnode'],
]:
    corrs = usedata.groupby(['sample', 'fiber_diam', 'contact', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
    corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    corrs['contact'] = pd.Categorical(corrs['contact'], categories=['3D', 'cathodic', 'center', 'anodic'], ordered=True)
    plt.figure()
    # g = sns.FacetGrid(data=corrs,col='contact')
    # g.map_dataframe(sns.stripplot, hue='Sample', y='correlation',dodge=True, palette='colorblind')
    # g.map_dataframe(sns.boxplot,y='correlation', boxprops={'facecolor':'None'},whis=100)
    sns.swarmplot(data=corrs, palette='RdPu', y='contact', hue='fiber_diam', x='correlation', dodge=True, s=6)
    sns.boxplot(data=corrs, y='contact', x='correlation', boxprops={'facecolor': 'None'}, whis=100)
    # plt.subplots_adjust(top=0.8)
    plt.title(f'Correlation between \n{comparison[0]} and {comparison[1]}', pad=25)
    plt.gca().set_xlim([-1, 1])
    plt.legend(title="D (μm)", bbox_to_anchor=(1, 1))
    plt.ylabel('')
    # plt.gca().set_yticklabels('')
    plt.gcf().set_size_inches([6, 5])
    plt.yticks(
        [0, 1, 2, 3], ['true-3D', 'cathodic-leading\nextrusion', 'center\nextrusion', 'anodic-leading\nextrusion']
    )
    # plt.ylabel('extrusion slice')
    # plt.title('')

    # # make lmplot to accompany
    # sns.lmplot(
    #     data=usedata.query('fiber_diam in [3,13]'),
    #     x=comparison[0],
    #     y=comparison[1],
    #     hue='nerve_label',
    #     col='contact',
    #     col_order=['cathodic', 'center', 'anodic'],
    #     palette='colorblind',
    #     row='fiber_diam',
    #     row_order=[3, 13],
    #     facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
    #     scatter_kws={'linewidth': 1, 'edgecolor': 'k'},
    # )
# %% imthera moreslice
imdata = pd.read_csv("thresh_unmatched_sim17_immynew.csv")
imdata = addpwfd(imdata, '3')
imdata['slice'] = (
    imdata['sample'].astype(str).str[2].replace({'2': 'cathodic', '0': 'anodic', '1': 'center'}).drop(columns='contact')
)
diams = {0: 3, 1: 13}
imdata['fiber_diam'] = imdata['fiberset_index'].replace(diams)

# inners need to match the center slice
innerid = {}
for inner in pd.unique(imdata.inner):
    innerid[inner] = pd.unique(imdata.query(f'inner == {inner} and type=="2D"')['master_fiber_index'])
    # get all rows where 3D and master fiber index is in the innerid
    imdata.loc[(imdata['type'] == '3D') & (imdata['master_fiber_index'].isin(innerid[inner])), 'inner'] = inner
imdata['type'] = imdata['type'].replace({'2D': 'extrusion', '3D': 'true-3D'})
im_matched = datamatch_merge(
    imdata.query('type=="extrusion"'),
    imdata.query('type=="true-3D"'),
    'threshold',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='type')


# Define a function to map the values
def categorize_contact(value):
    if value in [2, 5]:
        return '-2mm'
    elif value in [1, 4]:
        return 'center'
    elif value in [0, 3]:
        return '+2mm'
    else:
        raise ValueError('oops')


# Apply the function to create a new column
im_matched['srcpos'] = im_matched['active_src_index'].apply(categorize_contact)
im_matched['srcpos'] = pd.Categorical(im_matched['srcpos'], categories=['-2mm', 'center', '+2mm'], ordered=True)
# %% imthera unity
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
# usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
g = sns.relplot(
    data=im_matched.rename(columns={'nerve_label': 'Sample'}).query('fiber_diam ==3'),
    kind='scatter',
    col='active_src_index',
    row='slice',
    x='threshold',
    y='threshold3d',
    # hue='Sample',
    color='white',
    # s=20,
    palette='rainbow',
    facet_kws={'sharex': False, 'sharey': 'row', 'margin_titles': True},
    edgecolor='black',
    linewidth=1,
    alpha=1,
    hue='inner',
    legend=False,
)

for ax in g.axes.ravel():
    lim = np.amax([ax.get_xlim()[1], ax.get_ylim()[1]])
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='unity line')
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])
    ax.set_aspect('equal')
rownames(g, row_template='true-3D Threshold (mA)\nD = {row_name} μm')
plt.legend(loc='lower right')
g.set_xlabels('extrusion Threshold (mA)')
g.set_titles(col_template='Active contact: {col_name}', row_template='')
plt.subplots_adjust(wspace=-0.6)
# TODO need to edit all following analyses to support new slices
elceeceecee = []

for diam, contact in zip([3] * 6 + [13] * 6, [0, 1, 2, 3, 4, 5] * 2):
    for sli in im_matched.slice.unique():
        rdata = im_matched.query(f'fiber_diam==@diam and active_src_index=={contact} and slice==@sli')
        assert len(rdata) > 0
        r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
        print(f'{contact} {diam} μm: {r ** 2:.2f}')
        shortdat = rdata
        data2d = shortdat.sort_values('threshold').master_fiber_index
        data3d = shortdat.sort_values('threshold3d').master_fiber_index
        rc = compute_reorder_cost(list(data2d), list(data3d))
        elceeceecee.append(
            dict(diam=diam, contact=contact, value=r**2, slice=sli, ar=rc, srcpos=rdata.srcpos.unique()[0])
        )
elcdata = pd.DataFrame(elceeceecee)
elcdata['srcpos'] = pd.Categorical(elcdata['srcpos'], categories=['-2mm', 'center', '+2mm'], ordered=True)
elcdata['slice'].replace({"cathodic": "+4mm", "anodic": "-4mm"}, inplace=True)
plt.figure()
g = sns.lineplot(data=elcdata, style='diam', y='value', x='contact', hue='slice')
print(elcdata.groupby(['diam']).mean().round(2))
print(elcdata.groupby(['diam']).std().round(2))
sns.move_legend(g, [1, 0])
plt.ylabel('CCC')
plt.xticks(range(6))
g = sns.catplot(
    data=elcdata, kind='bar', hue='diam', y='value', col='contact', x='slice', col_wrap=2, col_order=[0, 3, 1, 4, 2, 5]
)
g.set_ylabels('CCC')
g = sns.catplot(
    data=elcdata, kind='bar', hue='diam', y='ar', col='contact', x='slice', col_wrap=2, col_order=[0, 3, 1, 4, 2, 5]
)
g.set_ylabels('AR')
g = sns.relplot(data=elcdata, kind='line', y='value', x='srcpos', hue='slice', errorbar=('se', 1))
plt.ylabel('$CCC^2$ ($R^2$)')
plt.xlabel('Source Position')
g.legend.set_title('Slice\nposition')
g = sns.relplot(data=elcdata, kind='line', y='ar', x='srcpos', hue='slice', errorbar=('se', 1))
plt.ylabel('Activation Reordering')
plt.xlabel('Source Position')
g.legend.set_title('Slice\nposition')


# %% imthera contact config
sns.set(font_scale=1.5, style='whitegrid')
imdata['active_src_index'] = pd.Categorical(imdata['active_src_index'], categories=[0, 3, 1, 4, 2, 5], ordered=True)

g = sns.catplot(
    data=imdata.query('fiber_diam in [3]'),
    x='type',
    y='threshold',
    kind='strip',
    # row='fiber_diam',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    hue='inner',
    col='active_src_index',
    linewidth=0.75,
    s=30,
    dodge=True,
    col_wrap=2,
    palette='rainbow',
    edgecolor='k',
)
g.set_ylabels('Threshold (mA)')
g.set_xlabels('')
g.set_titles(col_template='Active contact: {col_name}', row_template='D: {row_name} μm')
g.fig.set_size_inches(10, 10)
# %% imthera contact config lineplot
sns.set(font_scale=1.5, style='whitegrid')
imdata['active_src_index'] = pd.Categorical(imdata['active_src_index'], categories=[0, 3, 1, 4, 2, 5], ordered=True)

g = sns.relplot(
    data=imdata.query('fiber_diam in [3]'),
    x='inner',
    y='threshold',
    kind='line',
    # row='fiber_diam',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    hue='type',
    col='active_src_index',
    linewidth=0.75,
    # s=30,
    # dodge=True,
    col_wrap=2,
    palette=pal2d3d,
    # edgecolor='k',
    errorbar=('pi', 100),
)
g.set_ylabels('Threshold (mA)')
g.set_xlabels('fascicle')
g.legend.set_title('')
g.set_titles(col_template='Active contact: {col_name}', row_template='D: {row_name} μm')
g.fig.set_size_inches(10, 10)
# %% activation order imthera
# dose response
imdr = imdata.copy()
imdr['active_src_index'] = pd.Categorical(imdata['active_src_index'], categories=[0, 1, 2, 3, 4, 5], ordered=True)

imdr = calculate_dose_response(
    imdr,
    'threshold',
    'percent_activated',
    grouping_columns=['type', 'sample', 'fiber_diam', 'sim', 'active_src_index'],
)
imdr.sort_values('type', inplace=True)

# now do
newim = imdr.copy()
newim = newim.sort_values('type')
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
newim['percent_activated2d'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newim.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newim.query(
        f'type == "true-3D" and fiber_diam == @row.fiber_diam and master_fiber_index == @row.master_fiber_index and active_src_index==@row.active_src_index'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newim.loc[row.Index, 'percent_activated3d'] = val
newim['type'] = pd.Categorical(newim['type'], ordered=True, categories=['true-3D', 'extrusion'])
# plot
sns.set(font_scale=1.5, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    data=newim.query('fiber_diam==3'),
    y='percent_activated',
    x='type',
    units='nerve_label',
    col='active_src_index',
    palette='plasma',
    hue='percent_activated3d',
    # hue='inner',
    estimator=None,
    linewidth=0,
    facet_kws={'margin_titles': True},
    s=25,
)
plt.subplots_adjust(top=0.87)
# plt.suptitle(stringdat, x=0.37)
g.set_titles(row_template='', col_template='Active contact: {col_name}')
g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion fibers activated')
g.set_xlabels('')
g.legend.remove()
norm = plt.Normalize(0, 1)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])

# Remove the legend and add a colorbar
g.figure.colorbar(
    sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.8, label='Proportion extrusion\nfibers activated', pad=0.025
).ax.yaxis.set_ticks_position('left')

g.fig.set_size_inches(18, 6)
# %% reorder cost for imthera
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

scores = []
for cc in range(6):
    shortdat = imdata.query(f'fiber_diam==3 and active_src_index=={cc}')
    data2d = shortdat.query('type=="extrusion"').sort_values('threshold').master_fiber_index
    data3d = shortdat.query('type=="true-3D"').sort_values('threshold').master_fiber_index
    rc = compute_reorder_cost(list(data2d), list(data3d))
    scores.append({'active_src_index': cc, 'score2d3d': rc})
plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.barplot(data=scoredat, y='score2d3d', x='active_src_index', color='k')
ax.set_xlabel('Active Contact')
plt.ylabel('Activation Reordering')
plt.ylim(0, 0.5)
# %% reorder cost for imthera
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

scores = []
for cc in range(6):
    shortdat = imdata.query(f'fiber_diam==3 and active_src_index=={cc}')
    data2d = shortdat.query('type=="extrusion"')['threshold']
    data3d = shortdat.query('type=="true-3D"')['threshold']
    ccc = concordance_correlation_coefficient(data3d, data2d)
    scores.append({'active_src_index': cc, 'score2d3d': ccc})
plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.barplot(data=scoredat, y='score2d3d', x='active_src_index', color='k')
ax.set_xlabel('Active Contact')
plt.ylabel('CCC')

# %% reorder cost for imthera fascicle specific
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

mean_im = imdata.groupby(['active_src_index', 'inner', 'type', 'fiber_diam']).mean().reset_index()

scores = []
for cc in range(6):
    shortdat = mean_im.query(f'fiber_diam==3 and active_src_index=={cc}')
    data2d = shortdat.query('type=="extrusion"').sort_values('threshold').inner
    data3d = shortdat.query('type=="true-3D"').sort_values('threshold').inner
    rc = compute_reorder_cost(list(data2d), list(data3d))
    scores.append({'active_src_index': cc, 'score2d3d': rc})
plt.figure()
ax = sns.barplot(data=scoredat, y='score2d3d', x='active_src_index', color='k')
ax.set_xlabel('Active Contact')
plt.ylabel('Activation Reordering')
plt.ylim(0, 0.5)

# %% dose-response imthera
sns.set_style('white')
imdr = imdata.copy().rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
imdr = calculate_dose_response(
    imdr,
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim', 'active_src_index'],
)
imdr.sort_values('modeltype', inplace=True)
imdrmatch = datamatch_merge(
    imdr.query('modeltype=="extrusion"'),
    imdr.query('modeltype=="true-3D"'),
    'percent_activated',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='modeltype')

fig, axs = plt.subplots(1, 6, sharey=True, sharex=True)
for acsrc, ax in zip(sorted(pd.unique(imdr.active_src_index)), axs):
    # plt.figure()
    sns.set(font_scale=1.5, style='white')
    g = sns.scatterplot(
        data=imdr.query(f"fiber_diam in [3] and active_src_index==@acsrc"),
        y='percent_activated',
        x='threshold',
        hue='inner',
        palette='rainbow',
        linewidth=0,
        ax=ax,
        # legend=False
    )
    sns.lineplot(
        data=imdr.query(f"fiber_diam in [3] and active_src_index==@acsrc"),
        y='percent_activated',
        x='threshold',
        color='k',
        linewidth=2,
        estimator=None,
        units='modeltype',
        legend=False,
        zorder=1,
        ax=ax,
        style='modeltype',
    )

    ax.set_xlabel('Threshold (mA)')
    ax.set_ylim([0, 1])
    # plt.gcf().set_size_inches(8,4)
    # create legend, circle = extrusion, X= true-3D
    # create handles manually
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], linestyle='-', color='k', label='extrusion', linewidth=2),
        Line2D([0], [0], linestyle='--', color='k', label='true-3D', linewidth=2),
    ]
    legend_labels = ['extrusion', 'true-3D']
    g.legend(handles=legend_elements, labels=legend_labels, loc='lower right')
    ax.set_title(f'Active contact: {acsrc}')
fig.set_size_inches(24, 4)
axs[0].set_ylabel('Proportion Activated')
plt.subplots_adjust(wspace=0)
# %% for each contact configuration, and fiber diameter, calculate the fascicle selectivity ratio for each inner
imdata.reset_index(inplace=True, drop=True)


def recruitment_cost(data):
    """Calculate recruitment cost for each inner.

    Recruitment cost is defined as the ratio of number of stimulated off-target fibers to total number of off-target fibers.
    From https://iopscience.iop.org/article/10.1088/1741-2560/10/3/036010
    """
    for inner in pd.unique(data["inner"]):
        # get threshold for inner
        inner_data = data.query(f"inner == {inner}")['threshold']
        assert len(data) > 200 & len(data) < 250
        inner_thresh = np.amax(inner_data)
        # get all off-target fiber thresholds
        off_thresh = data.query(f"inner != '{inner}'")['threshold']
        # calculate recruitment cost
        cost = np.sum(off_thresh <= inner_thresh) / len(off_thresh)
        data.loc[data['inner'] == inner, 'RC'] = cost
        # fascicle selectivity ratio is 1-RC
        data.loc[data['inner'] == inner, 'FASR'] = 1 - cost
        fasr_dict = {
            'active_src_index': data['active_src_index'].iloc[0],
            'fiber_diam': data['fiber_diam'].iloc[0],
            'type': data['type'].iloc[0],
            'inner': inner,
            'RC': cost,
            'FASR': 1 - cost,
        }
        yield fasr_dict


imdatfasr = []
for contact_config in pd.unique(imdata['active_src_index']):
    for fiber_diam in pd.unique(imdata['fiber_diam']):
        for t in pd.unique(imdata['type']):
            imdatain = imdata.query(
                f'active_src_index == {contact_config} and fiber_diam == {fiber_diam} and type == "{t}"'
            )
            imdatfasr.extend(recruitment_cost(imdatain))

imdatfasr = pd.DataFrame(imdatfasr)
# %% plot FASR
sns.set(font_scale=2, style='white')
g = sns.relplot(
    data=imdatfasr.query('fiber_diam in [3]'),
    x='active_src_index',
    y='FASR',
    kind='line',
    col='inner',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    hue='inner',
    style='type',
    palette='rainbow',
    linewidth=3,
    col_wrap=7,
)
plt.ylim(0, 1)
g.set_xlabels('Active contact')
g.set_titles(col_template='fascicle: {col_name}')
plt.xticks(range(7))
# g.fig.set_size_inches([8, 10])
g.legend.remove()
g.fig.legend(handles=g.legend.legendHandles[-2:], labels=['extrusion', 'true-3D'], bbox_to_anchor=[0.94, 0.23])
# %% plot FASR
sns.set(font_scale=2, style='whitegrid')
imdatfnew = imdatfasr.query('fiber_diam in [3]')
imdatfnew = datamatch_merge(
    imdatfnew.query('type=="extrusion"'),
    imdatfnew.query('type=="true-3D"'),
    'FASR',
    merge_cols=["active_src_index", "inner"],
).drop(columns='type')
imdatfnew['FASR-diff'] = imdatfnew['FASR3d'] - imdatfnew['FASR']
g = sns.boxplot(
    data=imdatfnew,
    x='active_src_index',
    y='FASR-diff',
    # sharey=False,
    # hue='inner',
    # palette='rainbow',
    boxprops=dict(facecolor='none'),
    linewidth=3,
    # legend=False,
    whis=(0, 100),
)
g = sns.lineplot(
    data=imdatfnew,
    x='active_src_index',
    y='FASR-diff',
    # sharey=False,
    hue='inner',
    palette='rainbow',
    linewidth=2,
    legend=False,
    alpha=0.5,
)
plt.ylim(-1, 1)
plt.xlabel('Active contact')
plt.xticks(range(6))
# %% plot FASR
sns.set(font_scale=1.75, style='white')
g = sns.relplot(
    data=imdatfasr.query('fiber_diam in [3]'),
    x='inner',
    y='FASR',
    kind='line',
    col='active_src_index',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    style='type',
    palette='rainbow',
    linewidth=3,
    col_wrap=3,
)
plt.ylim(0, 1)
# %% calculate minthresh and activation order for all
newdefdr = defdr.copy()
anode = newdefdr.query('contact=="anodic" and modeltype=="extrusion"')
cathode = newdefdr.query('contact=="cathodic" and modeltype=="extrusion"')
match_min = datamatch_merge(
    anode,
    cathode,
    'threshold',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index", "deformation"],
).drop(columns='contact')
match_min['threshold'] = np.minimum(match_min['threshold'], match_min['threshold3d'])
matchmindr = calculate_dose_response(
    match_min.drop(columns='percent_activated'),
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim'],
)
matchmindr['contact'] = 'minthresh'
newdefdr = pd.concat([newdefdr, matchmindr])
newdefdr = newdefdr.query("deformation in ['Structural','Undeformed']").sort_values('modeltype')
newdefdr['percent_activated2d'] = np.nan
newdefdr.reset_index(drop=True, inplace=True)
# go through every row and for each fiber find the 2D activation percent
for row in newdefdr.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        'modeltype == "true-3D" and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and master_fiber_index == @row.master_fiber_index and deformation==@row.deformation'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newdefdr.loc[row.Index, 'percent_activated3d'] = val
    newdefdr.loc[row.Index, 'threshold3d'] = thisdat.threshold.values[0]
newdefdr['contact'] = pd.Categorical(
    newdefdr['contact'], categories=["3D", "cathodic", "center", "anodic", "minthresh"], ordered=True
)
theccc = newdefdr.copy()
# %% calc CCC and AR
sns.set(font_scale=1, style='whitegrid')
datahere = theccc.query('modeltype=="extrusion"')
scores = []
for cont in ["cathodic", "center", "anodic", "minthresh"]:
    for deformation in datahere.deformation.unique():
        for nerve in pd.unique(threshdat['nerve_label']):
            for n in [3, 13]:
                shortdat = datahere.query(
                    f'nerve_label=="{nerve}" and fiber_diam=={n} and contact==@cont and deformation==@deformation'
                )
                data2d = shortdat.sort_values('threshold').master_fiber_index
                data3d = shortdat.sort_values('threshold3d').master_fiber_index
                rc = compute_reorder_cost(list(data2d), list(data3d))
                data2d = shortdat.threshold
                data3d = shortdat.threshold3d
                ccc = concordance_correlation_coefficient(list(data2d), list(data3d))
                scores.append(
                    {
                        'sample': nerve,
                        'fiber_diam': n,
                        'score2d3d': rc,
                        'deformation': deformation,
                        'slice': cont,
                        'CCC': ccc,
                    }
                )
scoredat = pd.DataFrame(scores)
scoredat['fiber_diam'] = pd.Categorical(scoredat['fiber_diam'].astype(int), categories=[3, 13], ordered=True)
# barplot
g = sns.FacetGrid(data=scoredat, margin_titles=True, col="fiber_diam")
g.map_dataframe(sns.barplot, y='CCC', x='slice', hue='deformation', errorbar='se', palette='binary')
g.map_dataframe(
    sns.stripplot, y='CCC', x='slice', hue='deformation', dodge=True, color='black', edgecolor='white', linewidth=0.5
)

# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('CCC')
g.set_titles(col_template='D: {col_name} μm')
g.set_ylabels('CCC')
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
# plt.ylim(0, 0.6)
print(scoredat.groupby(['deformation']).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
g.fig.set_size_inches(8, 4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel('extrusion slice')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[2:], labels=labs[2:], bbox_to_anchor=[0.5, 1.3], ncol=2)
rownames(g, '$CCC^2$ -  D={row_name} μm')
# %% now plot
sns.set(font_scale=1, style='whitegrid')
g = sns.FacetGrid(data=scoredat, margin_titles=True, col="fiber_diam")
g.map_dataframe(sns.scatterplot, y='score2d3d', x='slice', style='deformation', hue='sample', palette=defpal, size=10)
g.map_dataframe(sns.lineplot, y='score2d3d', x='slice', style='deformation', hue='sample', palette=defpal)
# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('Activation Reordering')
g.set_titles(col_template='D: {col_name} μm')
g.set_ylabels('Activation Reordering')
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
plt.ylim(0, 0.6)
print(scoredat.groupby(['deformation']).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
# g.fig.set_size_inches(8,4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel('extrusion slice')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[9:], labels=labs[9:], bbox_to_anchor=[1.2, 1])
# barplot
g = sns.FacetGrid(data=scoredat, margin_titles=True, col="fiber_diam")
g.map_dataframe(sns.barplot, y='score2d3d', x='slice', hue='deformation', errorbar='se', palette='binary')
g.map_dataframe(
    sns.stripplot,
    y='score2d3d',
    x='slice',
    hue='deformation',
    dodge=True,
    color='black',
    edgecolor='white',
    linewidth=0.5,
)

# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('Activation Reordering')
g.set_titles(col_template='D: {col_name} μm')
g.set_ylabels('Activation Reordering')
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
plt.ylim(0, 0.6)
print(scoredat.groupby(['deformation']).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
# g.fig.set_size_inches(8,4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel('extrusion slice')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[2:], labels=labs[2:], bbox_to_anchor=[0.5, 1.4], ncol=2)
# barplot
g = sns.FacetGrid(data=scoredat, margin_titles=True, col="fiber_diam")
g.map_dataframe(sns.barplot, y='CCC', x='slice', hue='deformation', errorbar='se', palette='binary')
g.map_dataframe(
    sns.stripplot, y='CCC', x='slice', hue='deformation', dodge=True, color='black', edgecolor='white', linewidth=0.5
)

# ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('CCC')
g.set_titles(col_template='D: {col_name} μm')
g.set_ylabels('CCC')
# g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
# plt.ylim(0, 0.6)
print(scoredat.groupby(['deformation']).median())
# plt.subplots_adjust(wspace=0)
# plt.xlim(-1,2)
# g.fig.set_size_inches(8,4)
for ax in g.axes.ravel():
    plt.sca(ax)
    ax.set_xlabel('extrusion slice')
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[2:], labels=labs[2:], bbox_to_anchor=[0.5, 1.4], ncol=2)
# %% MCT
imdata = pd.read_csv("thresh_unmatched_sim121_MCT.csv")
imdata = addpwfd(imdata, '121', infile='plotconfig_MCT')
imdata['contact'] = imdata['sample'].astype(str).str[2].replace({'2': 'cathodic'})

imdata['fiber_diam'] = imdata['fiber_diam'].astype(int)

# inners for the true-3D data need to match their extrusion counterpart
for row in imdata.query('type=="3D"').itertuples():
    inner = imdata.query(
        f'type=="2D" and fiber_diam=={row.fiber_diam} and master_fiber_index=={row.master_fiber_index} and active_src_index=={row.active_src_index} and nerve_label=="{row.nerve_label}"'
    )['inner']
    assert len(inner) == 1
    imdata.loc[row.Index, 'inner'] = inner.values[0]

imdata['type'] = imdata['type'].replace({'2D': 'extrusion', '3D': 'true-3D'})
im_matched = datamatch_merge(
    imdata.query('type=="extrusion"'),
    imdata.query('type=="true-3D"'),
    'threshold',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='type')
# %% MCT unity
sns.set(context='paper', font_scale=1, style='whitegrid')

g = sns.relplot(
    data=im_matched,
    kind='scatter',
    col='active_src_index',
    x='threshold',
    y='threshold3d',
    hue='nerve_label',
    row='fiber_diam',
    color='white',
    # s=20,
    palette=defpal,
    facet_kws={'sharex': False, 'sharey': 'row', 'margin_titles': True},
    edgecolor='black',
    linewidth=1,
    alpha=1,
    # hue='inner',
)

g.set_xlabels('extrusion Threshold (mA)')
g.set_titles(col_template='Active contact: {col_name}', row_template='D: {row_name} μm')
g.set_ylabels('true-3D Threshold (mA)')
for ax in g.axes.ravel():
    lim = np.amax([ax.get_xlim()[1], ax.get_ylim()[1]])
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='unity line')
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])
    ax.set_aspect('equal')
plt.legend(loc='lower right')
rownames(g, row_template='true-3D Threshold (mA)\nD = {row_name} μm')
g.fig.set_size_inches(10, 10)

elceeceecee = []
# get data for each facet and calculate concordance correlation coefficient
for data in g.facet_data():
    shortdat = data[1]
    data2d = shortdat['threshold']
    data3d = shortdat['threshold3d']
    ccc = concordance_correlation_coefficient(data3d, data2d)
    elceeceecee.append(dict(contact=data[0][1], diam=data[0][0], value=ccc**2))
elcdata = pd.DataFrame(elceeceecee)
# plot
plt.figure()
sns.lineplot(data=elcdata, style='diam', y='value', x='contact')

# now calculate again but separate by nerve label
elceeceecee = []
# get data for each facet and calculate concordance correlation coefficient
for data in g.facet_data():
    for nerve in pd.unique(data[1]['nerve_label']):
        shortdat = data[1].query(f'nerve_label=="{nerve}"')
        data2d = shortdat['threshold']
        data3d = shortdat['threshold3d']
        ccc = concordance_correlation_coefficient(data3d, data2d)
        elceeceecee.append(dict(contact=data[0][1], diam=data[0][0], value=ccc**2, sample=nerve))
elcdata = pd.DataFrame(elceeceecee)

# plot
plt.figure()
sns.lineplot(data=elcdata, style='diam', y='value', x='contact', hue='sample')
plt.ylabel('ccc')

# %% MCT activation order
# dose response
imdr = imdata.copy()
imdr['active_src_index'] = pd.Categorical(imdata['active_src_index'], categories=[0, 1, 2, 3], ordered=True)

imdr = calculate_dose_response(
    imdr,
    'threshold',
    'percent_activated',
    grouping_columns=['type', 'sample', 'fiber_diam', 'sim', 'active_src_index'],
)
imdr.sort_values('type', inplace=True)

# now do
newim = imdr.copy()
newim = newim.sort_values('type')
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
newim['percent_activated2d'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newim.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newim.query(
        f'type == "true-3D" and fiber_diam == @row.fiber_diam and master_fiber_index == @row.master_fiber_index and active_src_index==@row.active_src_index and nerve_label==@row.nerve_label'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newim.loc[row.Index, 'percent_activated3d'] = val
newim['type'] = pd.Categorical(newim['type'], ordered=True, categories=['true-3D', 'extrusion'])
# plot
sns.set(font_scale=1.5, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    data=newim.query('fiber_diam==3'),
    y='percent_activated',
    x='type',
    units='nerve_label',
    col='active_src_index',
    palette='plasma',
    hue='percent_activated3d',
    # hue='inner',
    estimator=None,
    linewidth=0,
    facet_kws={'margin_titles': True},
    s=25,
    row='nerve_label',
)
# reorder cost for MCT
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

scores = []
for cc in range(4):
    shortdat = imdata.query(f'fiber_diam==3 and active_src_index=={cc}')
    for nerve_label in pd.unique(shortdat['nerve_label']):
        finallythedat = shortdat.query(f'nerve_label=="{nerve_label}"')
        data2d = finallythedat.query('type=="extrusion"').sort_values('threshold').master_fiber_index
        data3d = finallythedat.query('type=="true-3D"').sort_values('threshold').master_fiber_index
        rc = compute_reorder_cost(list(data2d), list(data3d))
        scores.append({'active_src_index': cc, 'score2d3d': rc, 'nerve_label': nerve_label})
plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.lineplot(data=scoredat, y='score2d3d', x='active_src_index', hue='nerve_label', palette=defpal)
ax.set_xlabel('Active Contact')
plt.ylabel('Activation Reordering')

# %% dose response MCT
sns.set_style('white')
imdr = imdata.copy().rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
imdr = calculate_dose_response(
    imdr,
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim', 'active_src_index', 'nerve_label'],
)
imdr.sort_values('modeltype', inplace=True)
imdrmatch = datamatch_merge(
    imdr.query('modeltype=="extrusion"'),
    imdr.query('modeltype=="true-3D"'),
    'percent_activated',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='modeltype')

for nerve_label in pd.unique(imdrmatch['nerve_label']):
    thisplotdata = imdr.query(f'nerve_label=="{nerve_label}"')
    fig, axs = plt.subplots(1, 4, sharey=True, sharex=True)
    for acsrc, ax in zip(sorted(pd.unique(thisplotdata.active_src_index)), axs):
        # plt.figure()
        sns.set(font_scale=1.5, style='white')
        g = sns.scatterplot(
            data=thisplotdata.query(f"fiber_diam in [3] and active_src_index==@acsrc"),
            y='percent_activated',
            x='threshold',
            hue='inner',
            palette='rainbow',
            linewidth=0,
            ax=ax,
            # legend=False
        )
        sns.lineplot(
            data=thisplotdata.query(f"fiber_diam in [3] and active_src_index==@acsrc"),
            y='percent_activated',
            x='threshold',
            color='k',
            linewidth=2,
            estimator=None,
            units='modeltype',
            legend=False,
            zorder=1,
            ax=ax,
            style='modeltype',
        )

        ax.set_xlabel('Threshold (mA)')
        ax.set_ylim([0, 1])
        plt.gcf().set_size_inches(20, 4)
        # create legend, circle = extrusion, X= true-3D
        # create handles manually
        from matplotlib.lines import Line2D

        legend_elements = [
            Line2D([0], [0], linestyle='-', color='k', label='extrusion', linewidth=2),
            Line2D([0], [0], linestyle='--', color='k', label='true-3D', linewidth=2),
        ]
        legend_labels = ['extrusion', 'true-3D']
        g.legend(handles=legend_elements, labels=legend_labels, loc='lower right')
        ax.set_title(f'Active contact: {acsrc}')
    axs[0].set_ylabel('Proportion fibers activated\nNerve: ' + nerve_label)
    plt.subplots_adjust(wspace=0)

# %% calculate FASR
imdata.reset_index(inplace=True, drop=True)


def recruitment_cost(data):
    """Calculate recruitment cost for each inner.

    Recruitment cost is defined as the ratio of number of stimulated off-target fibers to total number of off-target fibers.
    From https://iopscience.iop.org/article/10.1088/1741-2560/10/3/036010
    """
    for inner in pd.unique(data["inner"]):
        # get threshold for inner
        inner_data = data.query(f"inner == {inner}")['threshold']
        assert len(data) > 200 & len(data) < 250
        inner_thresh = np.amax(inner_data)
        # get all off-target fiber thresholds
        off_thresh = data.query(f"inner != '{inner}'")['threshold']
        # calculate recruitment cost
        cost = np.sum(off_thresh <= inner_thresh) / len(off_thresh)
        data.loc[data['inner'] == inner, 'RC'] = cost
        # fascicle selectivity ratio is 1-RC
        data.loc[data['inner'] == inner, 'FASR'] = 1 - cost
        fasr_dict = {
            'active_src_index': data['active_src_index'].iloc[0],
            'fiber_diam': data['fiber_diam'].iloc[0],
            'type': data['type'].iloc[0],
            'inner': inner,
            'RC': cost,
            'FASR': 1 - cost,
            'nerve_label': data['nerve_label'].iloc[0],
        }
        yield fasr_dict


imdatfasr = []
for contact_config in pd.unique(imdata['active_src_index']):
    for fiber_diam in pd.unique(imdata['fiber_diam']):
        for t in pd.unique(imdata['type']):
            for nerve_label in pd.unique(imdata['nerve_label']):
                imdatain = imdata.query(
                    f'active_src_index == {contact_config} and fiber_diam == {fiber_diam} and type == "{t}" and nerve_label=="{nerve_label}"'
                )
                imdatfasr.extend(recruitment_cost(imdatain))

imdatfasr = pd.DataFrame(imdatfasr)
# %% plot FASR
sns.set(font_scale=2, style='white')
imdatfnew = imdatfasr.query('fiber_diam in [3]')
imdatfnew = datamatch_merge(
    imdatfnew.query('type=="extrusion"'),
    imdatfnew.query('type=="true-3D"'),
    'FASR',
    merge_cols=["active_src_index", "inner", "nerve_label"],
).drop(columns='type')
imdatfnew['FASR-diff'] = imdatfnew['FASR3d'] - imdatfnew['FASR']
for nerve_label in pd.unique(imdatfasr['nerve_label']):
    plt.figure()
    imdatplot = imdatfnew.query(f'nerve_label=="{nerve_label}"')
    g = sns.boxplot(
        data=imdatplot.query('fiber_diam in [3]'),
        x='active_src_index',
        y='FASR-diff',
        # sharey=False,
        # hue='inner',
        # palette='rainbow',
        boxprops=dict(facecolor='none'),
        linewidth=3,
        # legend=False,
        whis=(0, 100),
    )
    g = sns.lineplot(
        data=imdatplot.query('fiber_diam in [3]'),
        x='active_src_index',
        y='FASR-diff',
        # sharey=False,
        hue='inner',
        palette='rainbow',
        linewidth=2,
        legend=False,
        alpha=0.5,
    )
    plt.ylim(-1, 1)
    plt.xlabel('Active contact')
    plt.xticks(range(4))
    plt.title(f'Nerve: {nerve_label}')
# %% MCT swarmplot with fascicle as hue
sns.set(font_scale=1.5, style='whitegrid')
for nerve in pd.unique(imdata['nerve_label']):
    plt.figure()
    g = sns.catplot(
        kind='swarm',
        data=imdata.query(f'fiber_diam==3 and nerve_label=="{nerve}"'),
        y='threshold',
        x='type',
        units='nerve_label',
        col='active_src_index',
        palette='rainbow',
        hue='inner',
        # hue='inner',
        estimator=None,
        linewidth=0,
        facet_kws={'margin_titles': True},
        s=25,
    )
    plt.subplots_adjust(top=0.87)
    # plt.suptitle(stringdat, x=0.37)
    g.set_titles(row_template='', col_template='Active contact: {col_name}')
    g.axes[0][0].set_xlabel('')
    g.axes[0][0].set_ylabel(f'Threshold (mA)\nNerve:{nerve}')
    g.set_xlabels('')
    g.legend.remove()
    norm = plt.Normalize(0, 1)
    sm = plt.cm.ScalarMappable(cmap="rainbow", norm=norm)
    sm.set_array([])

    # Remove the legend and add a colorbar
    g.figure.colorbar(
        sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.8, label='Fascicle', pad=0.025
    ).ax.yaxis.set_ticks_position('left')

    g.fig.set_size_inches(18, 6)
