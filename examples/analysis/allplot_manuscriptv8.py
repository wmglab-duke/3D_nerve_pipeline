# -*- coding: utf-8 -*-
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
        data.loc[data['waveform_index'] == nsim, 'pulse_width'] = pulse_width
        data.loc[data['fiberset_index'] == nsim, 'fiber_diam'] = fiber_diam
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


def rownames(g, row_template):
    args = dict(row_var=g._row_var)
    # Draw the row titles on the left edge of the grid
    for i, row_name in enumerate(g.row_names):
        ax = g.axes[i, 0]
        args.update(dict(row_name=row_name))
        title = row_template.format(**args)
        ax.set_ylabel(title)


gogo = "initial"

cath_comparison = ["cathodic", '3D', 3]
center_comparison = ["center", '3D', 3]
an_comparison = ["anodic", '3D', 3]
main_comparison = center_comparison

simNUM = input('Gimme simNUM')
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
#%% inners need to match the cathodic leading contact
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
#%%
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
#%%
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
newdefdat = addpwfd(newdefdat, '3')
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
#%% dose-response
print('deformdr')
defdr = newdefdat.copy().rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
defdr = calculate_dose_response(
    defdr,
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim', 'deformation'],
)
sys.exit('prepdone')
#%%violinplot
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
panelA['panel'] = 'A'

panelC = defdr.query(f"deformation!='Undeformed' and contact in {main_comparison}")
panelC['panel'] = 'C'

panelB = defdr.query(f"modeltype=='3DM' and deformation!='3D-3D' and contact in {main_comparison}")
panelB['panel'] = 'B'

panelD = defdr.query(f"deformation!='Undeformed' and contact in {cath_comparison}")
panelD['panel'] = 'D'

alldr = pd.concat([panelA, panelB, panelC, panelD])
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
g.set_titles(row_template='', col_template='{col_name}')

g.set_ylabels('Proportion activated')
g.set_xlabels('Threshold (mA)')
# change the line width for the legend
for line in g.legend.get_lines():
    line.set_linewidth(4.0)
# %% threshold vals
for stringdat in ['Undeformed', '3D-3D']:
    # deformstate = deformation!="Undeformed"
    thisdat = datamatch_merge(
        newdefdat.query(f'type=="2DEM" and deformation=="{stringdat}"'),
        newdefdat.query(f'type=="3DM" and deformation=="{stringdat}"'),
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
            pe_res = pe(np.mean(a), np.mean(b))
            pemean.append({'level': level, 'fiber_diam': fiber_diam, 'pe': pe_res})

    pemean = pd.DataFrame(pemean)
    pemean["deformation"] = stringdat
    pemeans.append(pemean)

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
)
# replace labels
# g._legend.texts[0].set_text('')
# g._legend.texts[5].set_text('')
g.set_xlabels('Fiber Diameter (μm)')
g.set_titles(col_template='Deformation: {col_name}')
g.set_ylabels("Percent Difference")
plt.xticks([3, 5, 7, 9, 11, 13])
plt.ylim([0, None])
ylimpe = g.axes[0][0].get_ylim()
g.legend.set_title('')
sns.move_legend(g, 'lower center', bbox_to_anchor=(0.78, 0.6), ncol=1)
for line in g.legend.get_lines():
    line.set_linewidth(2)
for handle in g.legend.legend_handles:
    handle.set_markersize(10)
# allpemean = pd.concat(pemeans)
# g = sns.relplot(
#     data=allpemean, kind='line', col='deformation', style='level', y='pe', x='fiber_diam', markers=True, color='k'
# )
# g.legend.set_title('')
# g.set_xlabels('Fiber Diameter (μm)')
# g.set_titles(col_template='Deformation: {col_name}')
# g.set_ylabels("Percent Difference")
# plt.xticks([3, 5, 7, 9, 11, 13])
# plt.ylim(ylimpe)
# sns.move_legend(g, 'lower center', bbox_to_anchor=(0.75, 0.5), ncol=1)
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

# %% threshold variances and coefficient of variation intrafascicle and inter
sns.set(font_scale=1.75, style='whitegrid')
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
plt.xlabel('Fiber Diameter (μm)')
# plt.yscale('log')
plt.gca().get_legend().remove()
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
#%% recruitment cost:
sns.set(font_scale=1)
sns.set_style('whitegrid')
datahere = deftomatch.query('deformation in ["3D-3D","Undeformed"]')
for comp in [main_comparison, cath_comparison]:
    scores = []
    for deformation in datahere.deformation.unique():
        threshdat = datahere.query(f'contact in {comp} and deformation==@deformation')
        for nerve in pd.unique(threshdat['nerve_label']):
            for n in [3, 13]:
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
    g.map_dataframe(sns.boxplot, y='score2d3d', x='fiber_diam', boxprops={'facecolor': 'white'}, whis=float('Inf'))
    g.map_dataframe(sns.stripplot, y='score2d3d', x='fiber_diam', color='black', size=5, jitter=0.15)
    ax.set_xlabel('Fiber Diameter (μm)')
    plt.ylabel('AR')
    g.set_titles(col_template='{col_name}')
    g.set_ylabels('AR')
    g.set_xlabels('D (μm)')
    # plt.title(f'{deformation} - {comp[0]}')
    plt.ylim(0, 0.6)
# %% imthera
imdata = pd.read_csv("thresh_unmatched_sim17_immytoo.csv")
imdata = addpwfd(imdata, '3')
imdata['contact'] = imdata['sample'].astype(str).str[2].replace({'2': 'cathodic', '0': 'anodic', '1': 'center'})

# inners need to match the center slice
innerid = {}
for inner in pd.unique(imdata.inner):
    innerid[inner] = pd.unique(imdata.query(f'inner == {inner} and type=="2D"')['master_fiber_index'])
    # get all rows where 3D and master fiber index is in the innerid
    imdata.loc[(imdata['type'] == '3D') & (imdata['master_fiber_index'].isin(innerid[inner])), 'inner'] = inner
imdata['type'] = imdata['type'].replace({'2D': '2DEM', '3D': '3DM'})
# %% imthera contact config
sns.set(font_scale=1.75, style='whitegrid')
g = sns.catplot(
    data=imdata.query('fiber_diam in [3]'),
    x='type',
    y='threshold',
    kind='strip',
    # row='fiber_diam',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    hue='inner',
    palette='rainbow',
    col='active_src_index',
    linewidth=0.75,
    s=30,
    dodge=True,
    col_wrap=3,
)
g.set_ylabels('Threshold (mA)')
g.set_xlabels('')
g.set_titles(col_template='{col_name}', row_template='D: {row_name} μm')
g.fig.set_size_inches(10, 6)
# %% imthera contact config
sns.set(font_scale=1.75, style='whitegrid')
g = sns.catplot(
    data=imdata.query('fiber_diam in [13]'),
    x='type',
    y='threshold',
    kind='strip',
    # row='fiber_diam',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    hue='inner',
    palette='rainbow',
    col='active_src_index',
    linewidth=0.75,
    s=30,
    dodge=True,
    col_wrap=3,
)
g.set_ylabels('Threshold (mA)')
g.set_xlabels('')
g.set_titles(col_template='{col_name}', row_template='D: {row_name} μm')
g.fig.set_size_inches(10, 6)
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
sns.set(font_scale=1.75, style='white')
g = sns.relplot(
    data=imdatfasr.query('fiber_diam in [3,13]'),
    x='active_src_index',
    y='FASR',
    kind='line',
    row='fiber_diam',
    facet_kws=dict(margin_titles=True),
    # sharey=False,
    hue='inner',
    style='type',
    palette='rainbow',
    linewidth=3,
)
plt.ylim(0, 1)
g.set_xlabels('Contact Configuration')
g.set_titles(col_template='D: {col_name} μm')
plt.xticks(range(7))
g.fig.set_size_inches([8, 10])
g.set_titles(row_template='D: {row_name} μm')
