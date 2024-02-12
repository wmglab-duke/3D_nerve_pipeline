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


def addpwfd(data, sim):
    with open('examples/analysis/plotconfig.json') as f:
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
threshload['type'] = threshload['type'].replace({'2D': '2DEM', '3D': '3DM'})
threshload = addpwfd(threshload, str(simNUM))
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
# 3DM and 2DEM trhresholds as different col same row
matched = datamatch_merge(
    threshload.query('type=="2DEM"'),
    threshload.query('type=="3DM"'),
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
repeated['type'] = repeated['type'].replace({'threshold': '2DEM', 'threshold3d': '3DM'})
repeated['type'] = pd.Categorical(repeated['type'], categories=['2DEM', '3DM'], ordered=True)
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
    drdat.query('modeltype=="2DEM"'),
    drdat.query('modeltype=="3DM"'),
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
# where nerve label contains "def" deformation is "3D-3D"
threshes.loc[threshes['nerve_label'].str.contains('def'), 'deformation'] = '3D-3D'
# where nerve label contains "asc" deformation is "2D-3D"
threshes.loc[threshes['nerve_label'].str.contains('asc'), 'deformation'] = '2D-3D'
# else deformation is "none"
threshes.loc[threshes['deformation'].isna(), 'deformation'] = "Undeformed"
# strip "def" and "asc" from nerve labels
threshes['nerve_label'] = threshes['nerve_label'].str.replace('def', '').str.replace('asc', '')

newdefdat = threshes.copy()

# remove all nerve_label = '6L'
newdefdat = newdefdat[newdefdat['nerve_label'] != '6L']
newdefdat = newdefdat[newdefdat['nerve_label'] != '2R']

newdefdat['type'] = newdefdat['type'].replace({'2D': '2DEM', '3D': '3DM'})
newdefdat = addpwfd(newdefdat, simNUM)
newdefdat['fiber_diam'] = newdefdat['fiber_diam'].astype(int)
# contact is cathodic if the third digit of sample int is 2, anodic if 0
newdefdat['contact'] = (
    newdefdat['sample'].astype(str).str[2].replace({'2': 'cathodic', '1': 'center', '0': 'anodic', '3': '3D'})
)

# set deformation as ordered categorical
newdefdat['deformation'] = pd.Categorical(
    newdefdat['deformation'], categories=['Undeformed', '2D-3D', '3D-3D'], ordered=True
)
newdefdat['nerve_label'] = pd.Categorical(newdefdat['nerve_label'], categories=['2L', '3R', '5R', '6R'], ordered=True)
newdefdat['deformed'] = newdefdat['deformation'] != "Undeformed"
deftomatch = newdefdat.copy()

# remove unused colors from palette
defpal = [sns.color_palette('colorblind')[ind] for ind in [0, 2, 3, 5]]
defdefcomp = newdefdat.query('type=="3DM" and deformation != "2D-3D"')
defdefcomp['deformed'] = defdefcomp['deformation'] != "Undeformed"
defsamples = ["2L", "3R", "5R", "6R"]
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
newdef['deformation'] = pd.Categorical(newdef['deformation'], categories=['Undeformed', '2D-3D', '3D-3D'], ordered=True)
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
panelA['panel'] = 'Undeformed 2DEM vs. 3DM\n(cathodic-leading slice)'

panelB = defdr.query(f"modeltype=='3DM' and deformation!='3D-3D' and contact in {main_comparison}")
panelB['panel'] = '3DM Undeformed vs. Deformed'

panelC = defdr.query(f"deformation!='Undeformed' and contact in {main_comparison}")
panelC['panel'] = 'Deformed 2DEM vs. 3DM\n(center slice)'

panelD = defdr.query(f"deformation!='Undeformed' and contact in {cath_comparison}")
panelD['panel'] = 'Deformed 2DEM vs. 3DM\n(cathodic-leading slice)'

alldr = pd.concat(reversed([panelA, panelB, panelC, panelD]))
alldr['nerve_label'] = pd.Categorical(alldr['nerve_label'], categories=defsamples, ordered=True)
alldr['deformation'] = pd.Categorical(alldr['deformation'], categories=['Undeformed', '2D-3D', '3D-3D'], ordered=True)

sns.set(font_scale=1.75, style='whitegrid')

plt.figure()
g = sns.relplot(
    kind='line',
    col='nerve_label',
    style='deformation',
    data=alldr.query(f"fiber_diam in [3]"),
    y='percent_activated',
    x='threshold',
    units='nerve_label',
    hue='modeltype',
    palette=pal2d3d,
    estimator=None,
    linewidth=4,
    facet_kws={'sharex': False, 'margin_titles': True},
    row='panel',
)

g.legend.set_title('')

for ax in g.axes.ravel():
    ax.set_ylim([0, 1])
    ax.set_xlim([0, None])
rownames(g, row_template="Proportion Activated\n{row_name}", rotation=0, labelpad=150)
g.set_titles(row_template='', col_template="{col_name}")
plt.subplots_adjust(hspace=0.25)

g.set_xlabels('Threshold (mA)')
# change the line width for the legend
for line, l in zip(g.legend.get_lines(), g.legend.get_texts()):
    line.set_linewidth(4.0)
    if l.get_text() in ['deformation', 'modeltype']:
        l.set_text('')
# %% threshold vals
for stringdat in ['Undeformed', '3D-3D']:
    # deformstate = deformation!="Undeformed"
    thisdat = datamatch_merge(
        newdefdat.query(f'type=="2DEM" and deformation=="{stringdat}" and contact in {main_comparison}'),
        newdefdat.query(f'type=="3DM" and deformation=="{stringdat}" and contact in {main_comparison}'),
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
for stringdat in ['Undeformed', '2D-3D', '3D-3D']:
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
                    f"nerve_label == '{sam}' and level == '{level}' and type == '2DEM' and fiber_diam == {fiber_diam}"
                )['threshold'].values
                b = compiled_data.query(
                    f"nerve_label == '{sam}' and level == '{level}' and type == '3DM' and fiber_diam == {fiber_diam}"
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
            a = compiled_data.query(f"level == '{level}' and type == '2DEM' and fiber_diam == {fiber_diam}")[
                'threshold'
            ].values
            b = compiled_data.query(f"level == '{level}' and type == '3DM' and fiber_diam == {fiber_diam}")[
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

# %% change in threshold after deformation
diffs = []
# calculate the change in onset and saturation for each sample (compared between undeformed and 3D-3D deformed)
for nerve_label in pd.unique(onsets_sats["Undeformed"].nerve_label):
    for fiber_diam in pd.unique(onsets_sats["Undeformed"].fiber_diam):
        for level in ['onset', 'saturation']:
            undeformed = (
                onsets_sats["Undeformed"]
                .query(
                    f"type=='3DM' and nerve_label == '{nerve_label}' and fiber_diam == {fiber_diam} and level == '{level}'"
                )['threshold']
                .values[0]
            )
            deformed = (
                onsets_sats["3D-3D"]
                .query(
                    f"type=='3DM' and nerve_label == '{nerve_label}' and fiber_diam == {fiber_diam} and level == '{level}'"
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
# # subtract 1 from sample if type is 3DM
# threshdat['sample'] = threshdat.apply(lambda x: x['sample'] - 1 if x.type == "3DM" else x['sample'], axis=1)
data = threshload.query(f"fiber_diam in [3,13] and contact in {main_comparison} and nerve_label=='5R'")

sns.reset_orig()
sns.set(font_scale=1.25, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
# plt.suptitle('thresholds compared between 2DEM (cathodic) and 3DM')
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
# plt.suptitle('thresholds compared between 2DEM (cathodic) and 3DM')

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
# plt.suptitle('min-max normalized thresholds compared between 3 um and 13 um thresholds (cath 2DEMEM)')
g.set_titles(col_template='{col_name}', row_template='')
g.set_ylabels('Normalized Threshold')
# plt.savefig('matchsim.png', dpi=400)
# g.axes.ravel()[0].set_ylabel('Threshold (mA)\nModel: 2DEM')
# g.axes.ravel()[6].set_ylabel('Threshold (mA)\nModel: 3DM')
g.fig.set_size_inches(5, 5)
# %% organizatino compare nsim
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshdat = threshload.query(f'contact in {main_comparison}')

scores = []
for nerve in pd.unique(threshdat['nerve_label']):
    for n in [3, 13]:
        shortdat = threshdat.query(f'nerve_label=="{nerve}" and fiber_diam=={n}')
        data2d = shortdat.query('type=="2DEM"').sort_values('threshold').master_fiber_index
        data3d = shortdat.query('type=="3DM"').sort_values('threshold').master_fiber_index
        rc = compute_reorder_cost(list(data2d), list(data3d))
        scores.append({'sample': nerve, 'fiber_diam': n, 'score2d3d': rc})
plt.figure()
scoredat = pd.DataFrame(scores)
scoredat['fiber_diam'] = pd.Categorical(scoredat['fiber_diam'].astype(int), categories=[3, 13], ordered=True)
ax = sns.boxplot(data=scoredat, y='score2d3d', x='fiber_diam', boxprops={'facecolor': 'white'}, whis=10)
ax = sns.stripplot(data=scoredat, y='score2d3d', x='fiber_diam', color='black', size=5, jitter=0.25)
ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('AR')
plt.ylim(0, 0.5)
print(scoredat.groupby(['fiber_diam']).median())
# plt.gcf().set_size_inches([3, 5])
# %% organizatino compare type
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshdat = threshload.query(f'contact in {main_comparison}')

scores = []
for nerve in pd.unique(threshdat['nerve_label']):
    for model in ['2DEM', '3DM']:
        shortdat = threshdat.query(f'nerve_label=="{nerve}" and type=="{model}"')
        datasmol = shortdat.query('nsim==0').sort_values('threshold').master_fiber_index
        databeeg = shortdat.query('nsim==5').sort_values('threshold').master_fiber_index
        rc = compute_reorder_cost(list(datasmol), list(databeeg))
        scores.append({'sample': nerve, 'type': model, 'scoresmolbeeg': rc})
plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.boxplot(data=scoredat, y='scoresmolbeeg', x='type', boxprops={'facecolor': 'white'}, whis=10)
ax = sns.stripplot(data=scoredat, y='scoresmolbeeg', x='type', color='black', size=5, jitter=0.25)
plt.ylabel('AR')
plt.xlabel('')
plt.ylim(0, 0.5)
print(scoredat.groupby(['type']).median())
# plt.gcf().set_size_inches([3, 5])
# %% activation order
newdefdr = defdr.copy()
newdefdr = newdefdr.query("'5R' in nerve_label and deformation!='2D-3D'").sort_values('modeltype')
newdefdr['nerve_label'] = pd.Categorical(newdefdr['nerve_label'], categories=['5R'])
newdefdr['deformation'] = pd.Categorical(newdefdr['deformation'], categories=['Undeformed', '3D-3D'])
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
newdefdr['percent_activated2d'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for i, row in newdefdr.iterrows():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        f'modeltype == "2DEM" and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim and master_fiber_index == @row.master_fiber_index and contact in {main_comparison} and deformation==@row.deformation'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newdefdr.loc[i, 'percent_activated2d'] = val
# plot
sns.set(font_scale=1.5, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    row='fiber_diam',
    data=newdefdr.query(f"fiber_diam in [3] and contact in {main_comparison}"),
    y='percent_activated',
    x='modeltype',
    units='nerve_label',
    col='deformation',
    palette='plasma',
    hue='percent_activated2d',
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

# Remove the legend and add a colorbar
g.figure.colorbar(
    sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.8, label='Proportion 2DEM\nfibers activated', pad=0.1
).ax.yaxis.set_ticks_position('left')

# g._legend.set_title('')
# %% recruitment cost:
sns.set(font_scale=1)
sns.set_style('whitegrid')
datahere = deftomatch.query('deformation in ["3D-3D","Undeformed"]')
for comp in [main_comparison, cath_comparison]:
    scores = []
    for deformation in datahere.deformation.unique():
        threshdat = datahere.query(f'contact in {comp} and deformation==@deformation')
        for nerve in pd.unique(threshdat['nerve_label']):
            for n in [3]:
                shortdat = threshdat.query(f'nerve_label=="{nerve}" and fiber_diam=={n}')
                data2d = shortdat.query('type=="2DEM"').sort_values('threshold').master_fiber_index
                data3d = shortdat.query('type=="3DM"').sort_values('threshold').master_fiber_index
                rc = compute_reorder_cost(list(data2d), list(data3d))
                scores.append(
                    {'sample': nerve, 'fiber_diam': n, 'score2d3d': rc, 'deformation': deformation, 'slice': comp[0]}
                )
    scoredat = pd.DataFrame(scores)
    scoredat['fiber_diam'] = pd.Categorical(scoredat['fiber_diam'].astype(int), categories=[3, 13], ordered=True)
    g = sns.FacetGrid(data=scoredat, col='deformation', margin_titles=True)
    g.map_dataframe(sns.boxplot, y='score2d3d', boxprops={'facecolor': 'white'}, whis=float('Inf'))
    g.map_dataframe(sns.stripplot, y='score2d3d', color='black', size=5, jitter=0.15)
    # ax.set_xlabel('Fiber Diameter (μm)')
    plt.ylabel('AR')
    g.set_titles(col_template='{col_name}')
    g.set_ylabels('AR')
    # g.set_xlabels('D (μm)')
    # plt.title(f'{deformation} - {comp[0]}')
    plt.ylim(0, 0.5)
    print(scoredat.groupby(['deformation']).median())
# %% remake the above with individual calls of histplot
sns.set(font_scale=1.75, style='whitegrid')
newthreshz = threshload.copy()
newthreshz['activation_zpos'] = newthreshz['activation_zpos'] / 10000
fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
for nerve_label in pd.unique(newthreshz.nerve_label):
    for ax, modeltype in zip(axs, ["2DEM", "3DM"]):
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
# plot one one line out to max of 3DM

for diam, pos, ax in zip([3, 13], (0.2, 0.8), g.axes.ravel()):
    rdata = nsimdata.query(f'fiber_diam=={diam}')
    r, p = pearsonr(rdata.threshold, rdata.threshold3d)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} μm: {r ** 2:.2f}')
    ax.set_title(f'{diam} μm')
    ax.plot([0, rdata.threshold.max()], [0, rdata.threshold.max()], '--k', linewidth=2, label='1:1 line')
    ax.set_xlabel('2DEM Threshold (mA)')
    ax.set_yticks(ax.get_xticks())
    ax.set_aspect('equal', 'box')
    ax.set_xlim([0, None])
    ax.set_ylim([0, ax.get_xlim()[1]])
g.axes.ravel()[0].set_ylabel('3DM Threshold (mA)')
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
plt.ylabel('Threshold CoV')
# plt.xlabel('Fiber Diameter (μm)')
# plt.yscale('log')
# plt.gca().get_legend().remove()
g.get_legend().set_title('D (μm)')
plt.gcf().set_size_inches([6, 4])
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
imdata['type'] = imdata['type'].replace({'2D': '2DEM', '3D': '3DM'})
im_matched = datamatch_merge(
    imdata.query('type=="2DEM"'),
    imdata.query('type=="3DM"'),
    'threshold',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='type')
# %% imthera 1:1
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
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='1:1 line')
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])
    ax.set_aspect('equal')
rownames(g, row_template='3DM Threshold (mA)\nD = {row_name} μm')
plt.legend(loc='lower right')
g.set_xlabels('2DEM Threshold (mA)')
g.set_titles(col_template='Active contact: {col_name}', row_template='')
plt.subplots_adjust(wspace=-0.6)

elceeceecee = []

for diam, contact in zip([3] * 6 + [13] * 6, [0, 1, 2, 3, 4, 5] * 2):
    rdata = im_matched.query(f'fiber_diam==@diam and active_src_index=={contact}')
    assert len(rdata) > 0
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    print(f'{contact} {diam} μm: {r ** 2:.2f}')
    elceeceecee.append(dict(diam=diam, contact=contact, value=r**2))
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
for i, row in newim.iterrows():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newim.query(
        f'type == "2DEM" and fiber_diam == @row.fiber_diam and master_fiber_index == @row.master_fiber_index and active_src_index==@row.active_src_index'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newim.loc[i, 'percent_activated2d'] = val
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
    hue='percent_activated2d',
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
    sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.8, label='Proportion 2DEM\nfibers activated', pad=0.05
).ax.yaxis.set_ticks_position('left')

# g._legend.set_title('')
# %% reorder cost for imthera
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

scores = []
for cc in range(6):
    shortdat = imdata.query(f'fiber_diam==3 and active_src_index=={cc}')
    data2d = shortdat.query('type=="2DEM"').sort_values('threshold').master_fiber_index
    data3d = shortdat.query('type=="3DM"').sort_values('threshold').master_fiber_index
    rc = compute_reorder_cost(list(data2d), list(data3d))
    scores.append({'active_src_index': cc, 'score2d3d': rc})
plt.figure()
scoredat = pd.DataFrame(scores)
ax = sns.barplot(data=scoredat, y='score2d3d', x='active_src_index', color='k')
ax.set_xlabel('Active Contact')
plt.ylabel('AR')
plt.ylim(0, 0.5)
# %% reorder cost for imthera fascicle specific
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

mean_im = imdata.groupby(['active_src_index', 'inner', 'type', 'fiber_diam']).mean().reset_index()

scores = []
for cc in range(6):
    shortdat = mean_im.query(f'fiber_diam==3 and active_src_index=={cc}')
    data2d = shortdat.query('type=="2DEM"').sort_values('threshold').inner
    data3d = shortdat.query('type=="3DM"').sort_values('threshold').inner
    rc = compute_reorder_cost(list(data2d), list(data3d))
    scores.append({'active_src_index': cc, 'score2d3d': rc})
plt.figure()
ax = sns.barplot(data=scoredat, y='score2d3d', x='active_src_index', color='k')
ax.set_xlabel('Active Contact')
plt.ylabel('AR')
plt.ylim(0, 0.5)

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
g.fig.legend(handles=g.legend.legendHandles[-2:], labels=['2DEM', '3DM'], bbox_to_anchor=[0.94, 0.23])
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
newthreshz = threshload.copy()
newthreshz['activation_zpos'] = newthreshz['activation_zpos'] / 10000
fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
for nerve_label in pd.unique(newthreshz.nerve_label):
    for ax, modeltype in zip(axs, ["2DEM", "3DM"]):
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
axs[0].set_title("2DEM-100%")
axs[1].set_title("3DM-100%")

# %% Correlation2DEM3DM
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = matched.rename(columns={'threshold': 'threshold2DEM'})
usedata.threshold3d = usedata.threshold3d.astype(float)
for comparison in [
    # ['threshold2DEM', 'threshold3d'],
    # ['threshold2DEM', 'peri_thk'],
    # ['threshold3d', 'peri_thk'],
    # ['threshold2DEM', 'minimum_efib_distance'],
    ['threshold2DEM', 'peak_second_diff'],
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
    # ['threshold', 'peak_second_diff'],
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
# %% Correlation2DEM3DM deformed
threeddefmatch = deftomatch.query('deformation=="3D-3D"')
deffinalmatch = datamatch_merge(
    threeddefmatch.query('type=="2DEM"'),
    threeddefmatch.query('type=="3DM"'),
    'threshold',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='type')
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = deffinalmatch.rename(columns={'threshold': 'threshold2DEM'})
for comparison in [
    ['threshold2DEM', 'threshold3d'],
    # ['threshold2DEM', 'peri_thk'],
    # ['threshold3d', 'peri_thk'],
    # # ['threshold2DEM', 'minimum_efib_distance'],
    # ['threshold2DEM', 'peak_second_diff'],
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
usedata = newdefdat.query('contact=="3D" and deformation=="3D-3D"')
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
thisdef = deftomatch.query('deformation!="2D-3D" and fiber_diam in [3,13]')
thisdef['deformed'] = thisdef['deformation'] == "3D-3D"

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
thisdef = deftomatch.query('deformation!="2D-3D" and fiber_diam in [3,13]')
thisdef['deformed'] = thisdef['deformation'] == "3D-3D"

comparecols = [
    'peak_second_diff',
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
g = sns.FacetGrid(data=alldf, col='type', row='fiber_diam', margin_titles=True)
g.map_dataframe(sns.pointplot, dodge=0.3, y='log_threshold', x='deformed', hue='name', palette='Set2')
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


# %%
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
for diam, pos, ax, deformation in zip([3, 3], (0.2, 0.8, 0.2, 0.8), g.axes.ravel(), ["Undeformed", "3D-3D"]):
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    lim = np.amax([ax.get_xlim()[1], ax.get_ylim()[1]])
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} μm: {r ** 2:.2f}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='1:1 line')
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])

    # ax.set_yticks(ax.get_xticks())
g.set_titles('D: {col_name} μm')
g.set_xlabels('2DEM Threshold (mA)')
g.set_ylabels('3DM Threshold (mA)')
g.set_titles(row_template='', col_template='Deformation: {col_name}')
# %% LCCC comparison maincomp
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
    [3, 13, 3, 13], (0.2, 0.8, 0.2, 0.8), g.axes.ravel(), ["Undeformed", "Undeformed", "3D-3D", "3D-3D"]
):
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    lim = ax.get_xlim()[1]
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} μm: {r ** 2:.2f}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='1:1 line')
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])

    # ax.set_yticks(ax.get_xticks())
g.set_titles(col_template='D: {col_name} μm')
g.set_xlabels('2DEM Threshold (mA)')
g.set_ylabels('3DM Threshold (mA)')
# %% remake the above with individual calls of histplot
sns.set(font_scale=1.75, style='whitegrid')
newthreshz = threshload.copy()
newthreshz['activation_zpos_oneten'] = newthreshz['activation_zpos_oneten'] / 10000
fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
for nerve_label in pd.unique(newthreshz.nerve_label):
    for ax, modeltype in zip(axs, ["2DEM", "3DM"]):
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
axs[0].set_title("2DEM-150%")
axs[1].set_title("3DM-150%")
# %% remake the above with individual calls of histplot
sns.set(font_scale=1.75, style='whitegrid')
newthreshz = threshload.copy()
newthreshz['activation_zpos'] = newthreshz['activation_zpos'] / 10000
fig, axs = plt.subplots(1, 2, sharex=False, sharey=True)
for nerve_label in pd.unique(newthreshz.nerve_label):
    for ax, modeltype in zip(axs, ["2DEM", "3DM"]):
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
axs[0].legend(loc='center left', bbox_to_anchor=(2.2, 0.5))
axs[0].set_xlim(reversed(axs[0].get_xlim()))
axs[0].set_ylabel("Activation Location (cm)")
axs[0].set_title("2DEM-110%")
axs[1].set_title("3DM-110%")
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
for stringdat in ['Undeformed', '2D-3D', '3D-3D']:
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
                        f"nerve_label == '{sam}' and level == '{level}' and type == '2DEM' and fiber_diam == {fiber_diam} and pulse_width == {pulse_width}"
                    )['threshold'].values
                    b = compiled_data.query(
                        f"nerve_label == '{sam}' and level == '{level}' and type == '3DM' and fiber_diam == {fiber_diam} and pulse_width == {pulse_width}"
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
    #         a = compiled_data.query(f"level == '{level}' and type == '2DEM' and fiber_diam == {fiber_diam")[
    #             'threshold'
    #         ].values
    #         b = compiled_data.query(f"level == '{level}' and type == '3DM' and fiber_diam == {fiber_diam}")[
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
# %% dose-response
# generate a median line across all samples
# so for 1%, find the first value for each sample which exceeds 1% and take the median of those
# then plot the median line for each fiber diameter
plt.figure()
sns.set(font_scale=1.5, style='white')
sns.lineplot(
    data=drdat.query(f"fiber_diam in [3] and contact in {main_comparison} and nerve_label =='5R'"),
    y='percent_activated',
    x='threshold',
    color='k',
    linewidth=1,
    estimator=None,
    units='modeltype',
    legend=False,
    zorder=1,
)

g = sns.scatterplot(
    data=drdat.query(f"fiber_diam in [3] and contact in {main_comparison} and nerve_label =='5R'"),
    y='percent_activated',
    x='threshold',
    hue='inner',
    palette='rainbow',
    linewidth=0,
    style='modeltype',
    # legend=False
)

plt.xlabel('Threshold (mA)')
plt.ylim([0, 1])
# plt.gcf().set_size_inches(8,4)
plt.ylabel('Proportion Activated')
# create legend, circle = 2DEM, X= 3DM
# create handles manually
from matplotlib.lines import Line2D

legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='2DEM', markerfacecolor='k', markersize=10),
    Line2D([0], [0], marker='X', color='w', label='3DM', markerfacecolor='k', markersize=10),
]
legend_labels = ['2DEM', '3DM']
g.legend(handles=legend_elements, labels=legend_labels, loc='lower right')
# %% activation order fiberdiam
newdefdr = defdr.copy()
newdefdr = newdefdr.query("'5R' in nerve_label and deformation!='2D-3D'").sort_values('modeltype')
newdefdr['nerve_label'] = pd.Categorical(newdefdr['nerve_label'], categories=['5R'])
newdefdr['deformation'] = pd.Categorical(newdefdr['deformation'], categories=['Undeformed', '3D-3D'])
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
newdefdr['percent_activated3'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for i, row in newdefdr.iterrows():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        f'fiber_diam == 3 and nerve_label == @row.nerve_label and modeltype == @row.modeltype and sim == @row.sim and master_fiber_index == @row.master_fiber_index and contact in {main_comparison} and deformation==@row.deformation'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newdefdr.loc[i, 'percent_activated3'] = val
# plot
sns.set(font_scale=1.5, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    # row='deformation',
    data=newdefdr.query(f"fiber_diam in [3,13] and contact in {main_comparison} and deformation=='Undeformed'"),
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
    sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.8, label='Proportion 2DEM\nfibers activated', pad=0.1
).ax.yaxis.set_ticks_position('left')
# %% LCCC comparison maincomp
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
for diam, pos, ax, deformation in zip([3, 3], (0.2, 0.8, 0.2, 0.8), g.axes.ravel(), ["Undeformed", "3D-3D"]):
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    lim[deformation] = np.amax([ax.get_xlim()[1], ax.get_ylim()[1]])
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} μm: {r ** 2:.2f}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim[deformation]], [0, lim[deformation]], '--k', linewidth=2, label='1:1 line')
    plt.legend()
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim[deformation]])
    ax.set_ylim([0, lim[deformation]])

    # ax.set_yticks(ax.get_xticks())
g.set_titles('D: {col_name} μm')
g.set_xlabels('2DEM Threshold (mA)')
g.set_ylabels('3DM Threshold (mA)')
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
for diam, pos, ax, deformation in zip([3, 3], (0.2, 0.8, 0.2, 0.8), g.axes.ravel(), ["Undeformed", "3D-3D"]):
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    r = concordance_correlation_coefficient(rdata.threshold3d, rdata.threshold)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} μm: {r ** 2:.2f}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim[deformation]], [0, lim[deformation]], '--k', linewidth=2, label='1:1 line')
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim[deformation]])
    ax.set_ylim([0, lim[deformation]])
    plt.legend()
    # ax.set_yticks(ax.get_xticks())
g.set_titles('D: {col_name} μm')
g.set_xlabels('2DEM Threshold (mA)')
g.set_ylabels('3DM Threshold (mA)')
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
                    f"nerve_label == '{sam}' and level == '{level}' and type == '2DEM' and fiber_diam == {fiber_diam} and pulse_width == {pulse_width}"
                )['threshold'].values
                b = compiled_data.query(
                    f"nerve_label == '{sam}' and level == '{level}' and type == '3DM' and fiber_diam == {fiber_diam} and pulse_width == {pulse_width}"
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
#         a = compiled_data.query(f"level == '{level}' and type == '2DEM' and fiber_diam == {fiber_diam")[
#             'threshold'
#         ].values
#         b = compiled_data.query(f"level == '{level}' and type == '3DM' and fiber_diam == {fiber_diam}")[
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
tcopy = threshload.copy()
tcopy['zres'] = (threshload.activation_zpos - threshload.peak_second_z) / (threshload.fiber_diam * 100)
g = sns.stripplot(data=tcopy, y='zres', x='type', hue='fiber_diam', dodge=True)
g.legend().set_title('D (μm)')
sns.move_legend(g, (1.05, 0))
plt.ylabel('Activation residual')
plt.xlabel('')
