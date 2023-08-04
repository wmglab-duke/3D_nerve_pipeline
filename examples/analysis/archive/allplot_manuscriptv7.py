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
from src.core.plotter import datamatch, datamatchlist

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
# TEMPORARY: remove 3R
# threshload.query("nerve_label != '3R'", inplace=True)
# sys.exit()
# END TEMPORARY
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
merged_df = pd.merge(df_3d, df_2d, on=['nerve_label', 'master_fiber_index', 'nsim'], suffixes=('_3d', '_2d'))

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
matched = datamatch(threshload.query('type=="2DEM"'), threshload.query('type=="3DM"'), 'threshold').drop(columns='type')
# add new EMsample column to matched dataframe, which is the nerve label plus the first letter of the contact type capitalized
matched['EMsample'] = matched['nerve_label'] + matched['contact'].str[0].str.upper()

# for matched boxplots
repeated = matched.melt(
    id_vars=[x for x in matched.columns if 'threshold' not in x],
    value_vars=['threshold3d', 'threshold'],
    var_name='type',
    value_name='threshold',
)
repeated['type'] = repeated['type'].replace({'threshold': '2DEM', 'threshold3d': '3DM'})
repeated['type'] = pd.Categorical(repeated['type'], categories=['2DEM', '3DM'], ordered=True)

# TODO: reactivate
multimatched = datamatchlist(
    threshload.query('type=="2DEM"'), threshload.query('type=="3DM"'), ['threshold', 'activation_zpos']
).drop(columns='type')
# %%
print('dose-response')
# dose response
drdat = threshload.copy()  # change this to repeated?
drdat['percent_activated'] = 0
drdat = drdat.rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
for i, row in drdat.iterrows():
    thisdat = drdat.query(
        'modeltype == @row.modeltype and samplenum == @row.samplenum and fiber_diam == @row.fiber_diam and sim == @row.sim'
    )
    # percent is number of thresholds less than or equal to this threshold divided by total number of thresholds
    drdat.loc[i, 'percent_activated'] = len(thisdat.query('threshold <= @row.threshold')) / len(thisdat)
drdat.sort_values('modeltype', inplace=True)

drmatch = datamatch(drdat.query('modeltype=="2DEM"'), drdat.query('modeltype=="3DM"'), 'percent_activated').drop(
    columns='modeltype'
)

sys.exit('prepdone')
# %% threshold violinplot
# generate boxplot of 2DEM and 3DM thresholds
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
g = sns.catplot(
    kind='violin',
    col='fiber_diam',
    data=threshload.query(f'fiber_diam in [3,13] and contact in {main_comparison}'),
    x='type',
    y='threshold',
    sharey=False,
    palette=pal2d3d,
    margin_titles=True,
    split=False,
)
for ax in g.axes.ravel():
    ax.set_xlabel('')
print(g.fig.get_size_inches())
g.fig.set_size_inches(5, 5, forward=True)
plt.subplots_adjust(top=0.9, wspace=0.4)
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# plt.suptitle('2DEM vs 3DM Thresholds for all samples')
g.axes.ravel()[0].set_title('3 μm')
g.axes.ravel()[1].set_title('13 μm')
plt.show()
# %% tortuosity example
sns.set(font_scale=1.25, style='whitegrid')
import matplotlib.pyplot as plt
import numpy as np


def calculate_tortuosity(x, y):
    total_distance = np.sum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
    straight_distance = np.sqrt((x[-1] - x[0]) ** 2 + (y[-1] - y[0]) ** 2)
    tortuosity = total_distance / straight_distance
    return tortuosity


def plot_tortuosity(non_directness, ax):
    num_points = 100
    x = np.zeros(num_points)
    y = np.linspace(0, 10, num_points)

    # Adding non-directness to the line
    deviation = 0.2 * non_directness
    x += np.random.uniform(-deviation, deviation, size=num_points)

    tortuosity = calculate_tortuosity(x, y)

    ax.plot(x, y, 'k', linewidth=2)
    ax.set_xlabel('X (a.u.)')
    ax.set_title(f"Tortuosity: {tortuosity:.2f}")
    ax.set_xlim(-1, 1)  # Set the same x limit for all plots
    ax.set_aspect('equal')
    ax.grid(False)


# Generating plots for different levels of non-directness
non_directness_values = [0, 0.1, 0.2]
# Create subplots
fig, axs = plt.subplots(1, len(non_directness_values), figsize=(6, 6), sharey=True)
axs[0].set_ylabel('Z (a.u.)')
for i, non_directness in enumerate(non_directness_values):
    plot_tortuosity(non_directness, axs[i])
plt.subplots_adjust(wspace=0)
plt.tight_layout()
plt.show()
# %% tortuosity lineplot
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
nsimdata = threshload.query('fiber_diam in [3,13] and type=="3DM"')
corrs = (
    nsimdata.groupby(['sample', 'fiber_diam', 'nerve_label', 'type'])['threshold', 'tortuosity'].corr().iloc[0::2, -1]
)
corrs = corrs.reset_index().rename(columns={'threshold': 'correlation'})
corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
corrs['tortuosity'] = corrs['tortuosity'] ** 2
means = corrs.groupby(['type', 'fiber_diam']).agg(np.mean)
# plot one one line out to max of 3DM
g = sns.lmplot(
    data=nsimdata.rename(columns={"nerve_label": "Sample"}),
    y='threshold',
    x='tortuosity',
    hue='Sample',
    palette='colorblind',
    col='fiber_diam',
    scatter_kws={'linewidths': 1, 'edgecolor': 'k', 's': 50},
    facet_kws={'margin_titles': True, 'sharey': False},
)
for ax, r in zip(g.axes.ravel(), means['tortuosity']):
    # ax.text(0.5, 0.9, f'\nmean r={r:.2f}', transform=ax.transAxes)
    print(r"$R^2$ = " + f"{r:.2f}")
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# plt.xlabel('3DM Threshold (mA)')
g.set_titles(col_template='{col_name} μm')
g.fig.set_size_inches([10, 5])
plt.show()
# %% dose-response
# generate a median line across all samples
# so for 1%, find the first value for each sample which exceeds 1% and take the median of those
# then plot the median line for each fiber diameter
newlinedat = []
for fiber_diam in [3, 13]:
    for percent in np.arange(0, 101, 1):
        for modeltype in pd.unique(drdat['modeltype']):
            datd = {"fiber_diam": [], "threshold": [], "percent": [], "model_type": []}
            thesevals = []
            for sampname in pd.unique(drdat.nerve_label):
                thisdat = drdat.query(
                    f'nerve_label=="{sampname}" and fiber_diam=={fiber_diam} and contact in {center_comparison} and modeltype=="{modeltype}"'
                )
                thisdat = thisdat.sort_values(by='threshold')
                thisdat = thisdat.loc[thisdat.percent_activated >= percent / 100]
                thisdat = thisdat.reset_index(drop=True)
                thesevals.append(thisdat.threshold.iloc[0])
            assert len(thesevals) == len(pd.unique(drdat['nerve_label']))
            datd["threshold"] = np.median(thesevals)
            datd["fiber_diam"] = fiber_diam
            datd["percent_activated"] = percent
            datd["model_type"] = modeltype
            newlinedat.append(datd)
newlinedat = pd.DataFrame(newlinedat)

plt.figure()
g = sns.relplot(
    kind='line',
    row='fiber_diam',
    data=drdat.query(f"fiber_diam in [3,13] and contact in {main_comparison}"),
    y='percent_activated',
    x='threshold',
    units='nerve_label',
    hue='modeltype',
    palette=pal2d3d,
    alpha=0.3,
    estimator=None,
    linewidth=2,
    facet_kws={'sharex': False},
)
for ax, fiber_diam in zip(g.axes.ravel(), [3, 13]):
    for color, modeltype in zip(pal2d3d, ['2DEM', '3DM']):
        x = np.array(newlinedat.query(f'fiber_diam=={fiber_diam} and model_type=="{modeltype}"').threshold)
        y = (
            np.array(newlinedat.query(f'fiber_diam=={fiber_diam} and model_type=="{modeltype}"').percent_activated)
            / 100
        )
        ax.plot(x, y, linewidth=3, color=color, linestyle='-')

g.legend.set_title('')
# sns.move_legend(g,[0.5,0.5])
for ax in g.axes.ravel():
    ax.set_xlabel('Threshold (mA)')
    ax.set_ylim([0, 1])
    ax.set_xlim([0, None])
g.set_titles(row_template='')

g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion of fibers activated\nFiber diameter: 3 μm')
g.axes[1][0].set_ylabel('Proportion of fibers activated\nFiber diameter: 13 μm')
g.axes[0][0].axhline(0.9, linewidth=2, color='k', linestyle='-', label='saturation')
g.axes[0][0].axhline(0.1, linewidth=2, color='k', linestyle='--', label='onset')
g.axes[1][0].axhline(0.9, linewidth=2, color='k', linestyle='-', label='saturation')
g.axes[1][0].axhline(0.1, linewidth=2, color='k', linestyle='--', label='onset')
# add text at 0.9 and 0.1 for each plot
g.axes[0][0].text(1.0, 0.87, 'saturation', fontsize=18, transform=g.axes[0][0].transAxes)
g.axes[0][0].text(1.0, 0.1, 'onset', fontsize=18, transform=g.axes[0][0].transAxes)
g.axes[1][0].text(1.0, 0.87, 'saturation', fontsize=18, transform=g.axes[1][0].transAxes)
g.axes[1][0].text(1.0, 0.1, 'onset', fontsize=18, transform=g.axes[1][0].transAxes)
g.fig.set_size_inches([7, 10])
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
# %% threshold vals
for nsim in pd.unique(matched.nsim):
    dadata = matched.query(f'nsim=={nsim} and contact in {main_comparison}')
    perc = sum(dadata.threshold > dadata.threshold3d) / len(dadata.threshold)
    print(f'Percent higher thresh for {nsim}:', round(perc, 3))
    res = np.abs(dadata.threshold3d - dadata.threshold)
    print(f'Mean abs residual for {nsim}:', round(np.mean(res), 3), 'sem:', round(sem(res), 3))
    print(
        f'Max,min,mean residual for {nsim}:',
        round(np.max(res), 3),
        round(np.min(res), 3),
        round(np.mean(res), 3),
        '+-',
        round(sem(res), 3),
    )
    mean_simthresh = np.mean(np.concatenate([dadata.threshold, dadata.threshold3d]))
    print(f'Mean sim thresh for {nsim}:', round(mean_simthresh, 3))
    # max, min, mean, sem as a percentage of mean_simthresh
    print(
        f'Max,min,mean,sem as a percentage of mean_simthresh for {nsim}:',
        round(np.max(res) / mean_simthresh, 3),
        round(np.min(res) / mean_simthresh, 3),
        round(np.mean(res) / mean_simthresh, 3),
        '+-',
        round(sem(res) / mean_simthresh, 3),
    )
    res = dadata.threshold3d - dadata.threshold
    print(f'Difference between means for {nsim}:', round(np.mean(res), 3), 'sem:', round(sem(res), 3))
# %% threshold variances and coefficient of variation intrafascicle and inter
sns.set(font_scale=1.75, style='whitegrid')
vardat = repeated.query(f'contact in {main_comparison}')
grouped = vardat.groupby(['EMsample', 'fiber_diam', 'type', 'inner'])
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

# %% all dose-response
newdr = drdat.copy()
sns.set(font_scale=2)
sns.set_style('whitegrid')
newdr['percent_activated2d'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for i, row in newdr.iterrows():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdr.query(
        f'modeltype == "2DEM" and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim and master_fiber_index == @row.master_fiber_index and contact in {main_comparison}'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newdr.loc[i, 'percent_activated2d'] = val
newdr.sort_values(by='samplenum')
# %%
sns.set(font_scale=1.75, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    row='fiber_diam',
    data=newdr.query(f"fiber_diam in [3,13] and contact in {main_comparison} and nerve_label=='5R'"),
    y='percent_activated',
    x='modeltype',
    units='nerve_label',
    col='nerve_label',
    palette='plasma',
    hue='percent_activated2d',
    # hue='inner',
    estimator=None,
    linewidth=0,
    facet_kws={'margin_titles': True},
    legend=False,
)

g.set_titles(row_template='', col_template='{col_name}')
g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion of fibers activated\nFiber diameter: 3 μm')
g.axes[1][0].set_ylabel('Proportion of fibers activated\nFiber diameter: 13 μm')
g.set_xlabels('')
# g._legend.set_title('Percentage of\n2D fibers activated')
# Remove the legend and add a colorbar
norm = plt.Normalize(0, 1)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])
g.figure.colorbar(
    sm, ax=g.axes.ravel().tolist(), aspect=14, shrink=0.6, label='Proportion 2DEM fibers activated', pad=0.15
).ax.yaxis.set_ticks_position('left')
g.fig.set_size_inches([10, 10])

# %% organization compare type with recruitment order
# remove all samples which dont end with a 2
sns.reset_orig()
sns.set(font_scale=1.25, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
# threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]
# # subtract 1 from sample if type is 3DM
# threshdat['sample'] = threshdat.apply(lambda x: x['sample'] - 1 if x.type == "3DM" else x['sample'], axis=1)
data = newdr.query(f"fiber_diam in [3,13] and contact in {main_comparison} and nerve_label=='5R'")

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
    x='modeltype',
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
    x='modeltype',
    y='threshold',
    palette='rainbow',
    hue='inner',
)

plt.subplots_adjust(top=0.9)
g.fig.set_size_inches(5, 5)
g.set_titles(col_template='{col_name}', row_template='')
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
axs[0].set_title("2DEM-100%")
axs[1].set_title("3DM-100%")
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
axs[0].legend(loc='center left', bbox_to_anchor=(2.2, 0.5))
axs[0].set_xlim(reversed(axs[0].get_xlim()))
axs[0].set_ylabel("Activation Location (cm)")
axs[0].set_title("2DEM-110%")
axs[1].set_title("3DM-110%")
# %% Correlation2DEM3DM
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = matched.rename(columns={'threshold': 'threshold2DEM'})
for comparison in [
    ['threshold2DEM', 'threshold3d'],
    ['threshold2DEM', 'peri_thk'],
    ['threshold3d', 'peri_thk'],
    ['threshold2DEM', 'minimum_efib_distance'],
    ['threshold2DEM', 'peak_second_diff'],
    ['peak_second_diff', 'peri_thk'],
    ['peak_second_diff', 'minimum_efib_distance'],
    ['peak_second_z', 'activation_zpos'],
    ['peak_second_diff_node', 'apnode'],
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
    ['threshold', 'peri_thk_act_site'],
    ['threshold', 'smallest_thk_under_cuff'],
    # ['threshold','peri_thk_act_site'],
    # ['smallest_thk_under_cuff', 'peri_thk_act_site'],
    # # TODO try smallest thk over whole cuffspan as well as over both cuff spans
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
# %%
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 400
sns.reset_orig()
tortuosities = [1, 1.01, 1.05, 1.1]
total_distance = 2  # meters
actual_distance = [total_distance / t for t in tortuosities]
diams = {10: 55.6, 11.5: 63.7, 12.8: 71.2, 14: 78.4, 15: 85.5, 16: 91.6}
eff = {t: {key: value / t for key, value in diams.items()} for t in tortuosities}

df = pd.DataFrame({'tortuosity': tortuosities, 'distance': actual_distance})

effdat = (
    pd.DataFrame(eff)
    .melt(ignore_index=False)
    .reset_index()
    .rename(columns={'index': 'diam', 'variable': 'tortuosity', 'value': 'cv'})
)

sns.lineplot(data=effdat, x='diam', y='cv', hue='tortuosity')
plt.xlabel("Fiber Diameter (μm)")
plt.ylabel("Conduction Velocity (m/s)")
plt.gca().set_aspect(0.12)

# sns.set(font_scale=1.75)
# sns.set_style('whitegrid')
plt.figure()
sns.histplot(data=threshload.query('type=="3DM"'), x='tortuosity', element='step')
tort = np.median(threshload.query('type=="3DM"').tortuosity)
plt.gca().set_aspect(0.00012)
print(f'Median tortuosity: {tort}')
# plt.title('Histogram of tortuosity values')
# %%
corrs = (
    nsimdata.groupby(['sample', 'fiber_diam', 'nerve_label', 'type'])['threshold', 'tortuosity'].corr().iloc[0::2, -1]
)
corrs = corrs.reset_index().rename(columns={'threshold': 'correlation'})
corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
corrs['R2'] = corrs['tortuosity'] ** 2
means = corrs.query('fiber_diam in [3,13]').groupby(['type', 'fiber_diam']).agg(R2=('R2', np.mean)).dropna()
# plot one one line out to max of 3DM
g = sns.lmplot(
    data=nsimdata.rename(columns={"nerve_label": "Sample"}).query('fiber_diam in [3,13]'),
    y='threshold',
    x='tortuosity',
    hue='Sample',
    palette='colorblind',
    col='fiber_diam',
    scatter_kws={'linewidths': 1, 'edgecolor': 'k', 's': 50},
    facet_kws={'margin_titles': True, 'sharey': False},
)
for ax, r in zip(g.axes.ravel(), means['R2']):
    # ax.text(0.5, 0.9, f'\nmean r={r:.2f}', transform=ax.transAxes)
    print(r"$R^2$ = " + f"{r:.2f}")
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# plt.xlabel('3DM Threshold (mA)')
g.set_titles(col_template='{col_name} μm')
g.fig.set_size_inches([10, 5])
plt.show()
# %% efib_distance lineplot
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
nsimdata = threshload.query(f'contact in {main_comparison}')
corrs = (
    nsimdata.groupby(['sample', 'fiber_diam', 'nerve_label', 'type'])[['threshold', 'minimum_efib_distance']]
    .corr()
    .iloc[0::2, -1]
)
corrs = corrs.reset_index().rename(columns={'threshold': 'correlation'})
corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
# print(corrs)
for fiber_diam in pd.unique(corrs['fiber_diam']):
    for nerve in pd.unique(nsimdata['nerve_label']):
        for typ in pd.unique(nsimdata['type']):
            # calculate slope of linear regression of threshold vs minimum_efib_distance
            this_data = nsimdata.query(f'fiber_diam=={fiber_diam} and nerve_label=="{nerve}" and type =="{typ}"')
            this_data['threshold'] = this_data['threshold'] / np.amax(this_data['threshold'])
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                this_data['minimum_efib_distance'] / 1000, this_data['threshold']
            )
            # add slope to corrs where fiber_diam and nerve_label match
            corrs.loc[
                (corrs['fiber_diam'] == fiber_diam) & (corrs['nerve_label'] == nerve) & (corrs['type'] == typ), 'slope'
            ] = slope
plt.gcf().set_size_inches(4, 4)
sns.swarmplot(data=corrs, y='slope', x='type', hue='fiber_diam', palette='RdPu', dodge=True)
plt.legend(bbox_to_anchor=[1, 1], title='D (μm)')
# plt.yscale('log')
plt.ylabel(r'slope $(\frac{mA/mA}{mm})$')
plt.xlabel('')
corrs['R2'] = corrs['minimum_efib_distance'] ** 2
means = (
    corrs.query('fiber_diam in [3,13]').groupby(['type', 'fiber_diam'])["minimum_efib_distance"].agg(np.mean).dropna()
)
print(means)
# plot one one line out to max of 3DM
g = sns.lmplot(
    data=nsimdata.rename(columns={"nerve_label": "Sample"}).query('fiber_diam in [3,13]'),
    y='threshold',
    x='minimum_efib_distance',
    hue='Sample',
    palette='colorblind',
    col='fiber_diam',
    row='type',
    scatter_kws={'linewidths': 1, 'edgecolor': 'k', 's': 50},
    facet_kws={'margin_titles': True, 'sharey': False},
)
for ax, r in zip(g.axes.ravel(), means):
    # ax.text(0.5, 0.9, f'\nmean r={r:.2f}', transform=ax.transAxes)
    print(r"$R^2$ = " + f"{r:.2f}")
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# plt.xlabel('3DM Threshold (mA)')
g.set_titles(col_template='{col_name} μm')
g.set_xlabels("Minimum Electrode-Fiber\nDistance (μm)")
g.set_ylabels('Threshold (mA)')
g.fig.set_size_inches([10, 10])
plt.show()
# %% peri_thk lineplot
val = 'peri_thk'
# val = 'peri_thk'
xlabel = "Activation Site Perineurium\nThickness(μm)"
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
nsimdata = threshload.query(f'contact in {main_comparison}')
# nsimdata = newdefdat.query(f'contact in {main_comparison} and deformation=="3D-3D"')
corrs = nsimdata.groupby(['sample', 'fiber_diam', 'nerve_label', 'type'])[['threshold', val]].corr().iloc[0::2, -1]
corrs = corrs.reset_index().rename(columns={'threshold': 'correlation'})
corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
# print(corrs)
for fiber_diam in pd.unique(corrs['fiber_diam']):
    for nerve in pd.unique(nsimdata['nerve_label']):
        for typ in pd.unique(nsimdata['type']):
            # calculate slope of linear regression of threshold vs peri_thk_act_site
            this_data = nsimdata.query(f'fiber_diam=={fiber_diam} and nerve_label=="{nerve}" and type =="{typ}"')
            this_data['threshold'] = this_data['threshold'] / np.amax(this_data['threshold'])
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                this_data[val] / 1000, this_data['threshold']
            )
            # add slope to corrs where fiber_diam and nerve_label match
            corrs.loc[
                (corrs['fiber_diam'] == fiber_diam) & (corrs['nerve_label'] == nerve) & (corrs['type'] == typ), 'slope'
            ] = slope
plt.gcf().set_size_inches(4, 4)
sns.swarmplot(data=corrs, y='slope', x='type', hue='fiber_diam', palette='RdPu', dodge=True)
plt.legend(bbox_to_anchor=[1, 1], title='D (μm)')
# plt.yscale('log')
plt.ylabel(r'slope $(\frac{mA/mA}{mm})$')
plt.xlabel('')
corrs['R2'] = corrs[val] ** 2
means = corrs.query('fiber_diam in [3,13]').groupby(['type', 'fiber_diam'])[val].agg(np.mean).dropna()
print(means)
# plot one one line out to max of 3DM
g = sns.lmplot(
    data=nsimdata.rename(columns={"nerve_label": "Sample"}).query('fiber_diam in [3,13]'),
    y='threshold',
    x='peri_thk_act_site',
    hue='Sample',
    palette='colorblind',
    col='fiber_diam',
    row='type',
    scatter_kws={'linewidths': 1, 'edgecolor': 'k', 's': 50},
    facet_kws={'margin_titles': True, 'sharey': False},
)
for ax, r in zip(g.axes.ravel(), means):
    # ax.text(0.5, 0.9, f'\nmean r={r:.2f}', transform=ax.transAxes)
    print(r"$R^2$ = " + f"{r:.2f}")
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# plt.xlabel('3DM Threshold (mA)')
g.set_titles(col_template='{col_name} μm')
g.set_xlabels(xlabel)
g.set_ylabels('Threshold (mA)')
g.fig.set_size_inches([10, 10])
plt.show()
# %% second diff
val = 'peak_second_diff'
xlabel = "Second Difference Ve"
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
nsimdata = threshload.query(f'contact in {main_comparison}')
# nsimdata = newdefdat.query(f'contact in {main_comparison} and deformation=="3D-3D"')
corrs = nsimdata.groupby(['sample', 'fiber_diam', 'nerve_label', 'type'])[['threshold', val]].corr().iloc[0::2, -1]
corrs = corrs.reset_index().rename(columns={'threshold': 'correlation'})
corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
# print(corrs)
for fiber_diam in pd.unique(corrs['fiber_diam']):
    for nerve in pd.unique(nsimdata['nerve_label']):
        for typ in pd.unique(nsimdata['type']):
            # calculate slope of linear regression of threshold vs peri_thk_act_site
            this_data = nsimdata.query(f'fiber_diam=={fiber_diam} and nerve_label=="{nerve}" and type =="{typ}"')
            this_data['threshold'] = this_data['threshold'] / np.amax(this_data['threshold'])
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                this_data[val] / 1000, this_data['threshold']
            )
            # add slope to corrs where fiber_diam and nerve_label match
            corrs.loc[
                (corrs['fiber_diam'] == fiber_diam) & (corrs['nerve_label'] == nerve) & (corrs['type'] == typ), 'slope'
            ] = slope
plt.gcf().set_size_inches(4, 4)
sns.swarmplot(data=corrs, y='slope', x='type', hue='fiber_diam', palette='RdPu', dodge=True)
plt.legend(bbox_to_anchor=[1, 1], title='D (μm)')
# plt.yscale('log')
plt.ylabel(r'slope $(\frac{mA/mA}{mm})$')
plt.xlabel('')
corrs['R2'] = corrs[val] ** 2
means = corrs.query('fiber_diam in [3,13]').groupby(['type', 'fiber_diam'])[val].agg(np.mean).dropna()
print(means)
# plot one one line out to max of 3DM
g = sns.lmplot(
    data=nsimdata.rename(columns={"nerve_label": "Sample"}).query('fiber_diam in [3,13]'),
    y='threshold',
    x='peri_thk_act_site',
    hue='Sample',
    palette='colorblind',
    col='fiber_diam',
    row='type',
    scatter_kws={'linewidths': 1, 'edgecolor': 'k', 's': 50},
    facet_kws={'margin_titles': True, 'sharey': False},
)
for ax, r in zip(g.axes.ravel(), means):
    # ax.text(0.5, 0.9, f'\nmean r={r:.2f}', transform=ax.transAxes)
    print(r"$R^2$ = " + f"{r:.2f}")
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# plt.xlabel('3DM Threshold (mA)')
g.set_titles(col_template='{col_name} μm')
g.set_xlabels(xlabel)
g.set_ylabels('Threshold (mA)')
g.fig.set_size_inches([10, 10])
plt.show()
# %% Percent Error
sns.reset_orig()
sns.set_style('whitegrid')
sns.set(font_scale=1.75, style='whitegrid')
# mpl.rcParams['figure.dpi'] = 400
# mpl.rcParams['font.size'] = 14
# sns.set(font_scale=1.75)
# apply pe to all rows of dataframe matched, with threshold3d as the correct value and threshold as the estimated value
matched['pe'] = matched.apply(lambda row: pe(row['threshold3d'], row['threshold']), axis=1)
# calculate difference between activation_zpos and activation_zpos3d
# multimatched['zdiff'] = multimatched['activation_zpos'] - multimatched['activation_zpos3d']
# multimatched['zdiff_abs'] = multimatched['zdiff'].abs()
plt.figure()
sns.barplot(data=matched, x='nerve_label', y='pe', hue='fiber_diam', errorbar='se', palette="RdPu")
# plt.title('Threshold Percent Error by sample and fiber diameter')
plt.legend(title='D (μm)', bbox_to_anchor=[1, 1], ncols=2)
plt.xlabel('')
plt.ylabel('Percent Difference (%)')
plt.gcf().set_size_inches([6, 5])
# plt.gca().set_aspect(0.04)
# sns.move_legend(plt.gca(), "upper left", bbox_to_anchor=(1, 1))
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
# %% stimfigs
fig, axs = plt.subplots(2, 2)
axs[0][0].plot([0, 1, 1, 2, 2, 3, 3, 4], [0, 0, -1, -1, 1, 1, 0, 0], 'k', linewidth=3)
axs[1][0].plot([0, 1, 1, 2, 2, 3, 3, 4], [0, 0, 1, 1, -1, -1, 0, 0], 'k', linewidth=3)
axs[0][1].plot([0, 1, 1, 2, 2, 3, 3, 4], [0, 0, -1, -1, 0, 0, 0, 0], 'k', linewidth=3)
axs[0][1].set_ylim(axs[0][0].get_ylim())
axs[1][1].plot([0, 1, 1, 2, 2, 3, 3, 4], [0, 0, 0, 0, 0, 0, 0, 0], 'k', linewidth=3)
axs[1][1].set_ylim(axs[0][0].get_ylim())
for ax in axs.ravel():
    ax.axis('off')
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
ax = sns.boxplot(data=scoredat, y='score2d3d', x='fiber_diam', boxprops={'facecolor': 'white'})
ax = sns.stripplot(data=scoredat, y='score2d3d', x='fiber_diam', color='black', size=5, jitter=0.25)
ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('AR')
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
ax = sns.boxplot(data=scoredat, y='scoresmolbeeg', x='type', boxprops={'facecolor': 'white'})
ax = sns.stripplot(data=scoredat, y='scoresmolbeeg', x='type', color='black', size=5, jitter=0.25)
plt.ylabel('AR')
plt.xlabel('')
# plt.gcf().set_size_inches([3, 5])
# %% peri-thk with 2DEM3DM
sns.set(font_scale=2)
sns.set_style('whitegrid')
nsimdata = threshload.query(f'fiber_diam in [3,13] and contact in {main_comparison}')
corrs = nsimdata.groupby(['sample', 'fiber_diam', 'nerve_label', 'type'])['threshold', 'peri_thk'].corr().iloc[0::2, -1]
corrs = corrs.reset_index().rename(columns={'threshold': 'correlation'})
corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
means = corrs.groupby(['type', 'fiber_diam']).agg(np.mean)
# plot one one line out to max of 3DM
g = sns.lmplot(
    data=nsimdata.rename(columns={"nerve_label": "Sample"}),
    y='threshold',
    x='peri_thk',
    hue='Sample',
    palette='colorblind',
    col='type',
    scatter_kws={'linewidths': 1, 'edgecolor': 'k', 's': 50},
    row='fiber_diam',
    facet_kws={'margin_titles': True, 'sharey': False},
)
g.set_ylabels('Threshold (mA)')
g.set_xlabels('Perineurium Thickness (μm)')
for ax, r, letter in zip(
    g.axes.ravel(order='F'), means['peri_thk'], ['A', '', 'B', '']
):  # TODO: these are not ordering correctly
    nl = '\n'
    ax.text(0, 0.9, fr'{letter}{nl}$\bar{{r}}$={r:.2f}', transform=ax.transAxes)
# plt.xlabel('3DM Threshold (mA)')
g.set_titles(col_template='{col_name}', row_template='Fiber Diameter = {row_name}')
# g.fig.set_size_inches(10, 6)
plt.show()
# %% Dose-response info all diams
plt.figure()
levels = {
    'onset': 10,
    'saturation': 90,
}
grouped = threshload.groupby(['sample', 'fiber_diam', 'sim', 'type', 'nerve_label', 'model', 'nsim', 'contact'])
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
    id_vars=['sample', 'fiber_diam', 'sim', 'type', 'nerve_label', 'model', 'nsim', 'contact'],
    value_vars=['onset', 'saturation'],
    var_name='level',
    value_name='threshold',
)

# set up facetgrid with nsim as row and level as columns
compiled_data.reset_index(inplace=True)
# set fiber_diam to category
compiled_data.type = compiled_data.type.astype('category')
# remove all rows where modulus sample with 0 is 0
compiled_data = compiled_data.query(f'contact in {main_comparison}')
# add a units column with unique number for each combination of fiber_diam and level
compiled_data['units'] = compiled_data.groupby(['fiber_diam', 'level', 'nerve_label']).ngroup()
compiled_data['fiber_diam'] = compiled_data['fiber_diam'].astype(int)
sns.set(font_scale=1.15, style='whitegrid')
g = sns.FacetGrid(
    compiled_data.query("fiber_diam in [3,13]").rename(
        columns={'fiber_diam': "Fiber Diameter (μm)", 'threshold': "Threshold (mA)"}
    ),
    col="level",
    row="Fiber Diameter (μm)",
    sharey=False,
    margin_titles=True,
)
g.map_dataframe(sns.swarmplot, x='type', y='Threshold (mA)', color='black')
g.map_dataframe(sns.lineplot, x='type', y='Threshold (mA)', units='units', hue='nerve_label', estimator=None, color='k')
plt.subplots_adjust(top=0.9)
g.set_titles(col_template='{col_name}', row_template='')
for ax in g.axes.ravel():
    ax.set_xlabel('')
g.axes[0][0].set_ylabel('Threshold (mA)\n3 μm fibers')
g.axes[1][0].set_ylabel('Threshold (mA)\n13 μm fibers')
plt.gcf().set_size_inches([5, 7])
plt.legend(title='Sample', bbox_to_anchor=(1.8, 1.5))

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
# plot percent error for each sample DR
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
plt.figure()
sns.lineplot(
    data=pes,
    style='level',
    y='pe',
    x='fiber_diam',
    markers=True,
    estimator=None,
    units='nerve_label',
    hue='nerve_label',
)
plt.xlabel('Fiber Diameter (μm)')
plt.ylabel('Percent Difference')
# plt.title('Percent Error for Onset and Saturation')
errorlim = plt.gca().get_ylim()
plt.xticks(pd.unique(compiled_data['fiber_diam']))
# plt.gca().get_legend().remove()
# plt.gcf().set_size_inches([7,5])
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[1:], labels=labs[1:], bbox_to_anchor=[1, 1])


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

# plot percent error for population
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
plt.figure()
sns.lineplot(data=pemean, style='level', y='pe', x='fiber_diam', markers=True, color='k')
plt.xlabel('Fiber Diameter (μm)')
plt.ylabel('Percent Difference')
# plt.title('Percent Error for pop Onset and Saturation')
plt.xticks(pd.unique(compiled_data['fiber_diam']))
plt.ylim(errorlim)
plt.legend(title='', ncols=2, bbox_to_anchor=[0.9, 1.2])
# %% BEGIN DEFORMATION ANALYSIS
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
deftomatch = newdefdat.copy()

# remove all aodic
newdefdat = newdefdat.query(f'contact in {main_comparison}')
# remove unused colors from palette
defpal = [sns.color_palette('colorblind')[ind] for ind in [0, 2, 3, 5]]
defdefcomp = newdefdat.query('type=="3DM" and deformation != "2D-3D"')
defdefcomp['deformed'] = defdefcomp['deformation'] != "Undeformed"


# %% dose-response
def calculate_dose_response(df, threshold_column, outcolumn, grouping_columns):
    def calculate_dose_response_per_sample(group):
        threshold_values = group[threshold_column].values
        activated_fibers = group[threshold_column].values[:, None] <= threshold_values
        dose_response = activated_fibers.sum(axis=0) / len(group)
        group[outcolumn] = dose_response
        return group

    return df.groupby(grouping_columns).apply(calculate_dose_response_per_sample)


defdr = newdefdat.copy().rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
defdr = calculate_dose_response(
    defdr,
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim', 'deformation'],
)
# %% tortuosities
sns.set(font_scale=1.75)
sns.set_style('white')
plt.figure()
sns.histplot(data=threshload.query('type=="3DM"'), x='cuff_tortuosity', element='step')
xl = plt.xlim()
plt.figure()
sns.histplot(data=newdefdat.query('type=="3DM" and deformation =="3D-3D"'), x='cuff_tortuosity', element='step')
plt.xlim(xl)

for nerve_label in pd.unique(threshload.nerve_label):
    print(nerve_label, np.median(threshload.query('nerve_label==@nerve_label and type=="3DM"')['cuff_tortuosity']))
for nerve_label in pd.unique(newdefdat.query('deformation=="Undeformed"').nerve_label):
    print(
        nerve_label,
        np.median(
            newdefdat.query('nerve_label==@nerve_label and type=="3DM" and deformation=="Undeformed"')[
                'cuff_tortuosity'
            ]
        ),
    )
for nerve_label in pd.unique(newdefdat.query('deformation=="3D-3D"').nerve_label):
    print(
        nerve_label,
        np.median(
            newdefdat.query('nerve_label==@nerve_label and type=="3DM" and deformation=="3D-3D"')['cuff_tortuosity']
        ),
    )
tort = np.median(newdefdat.query('type=="3DM" and deformation=="Undeformed"').cuff_tortuosity)
print(f'Median tortuosity: {tort}')
tort = np.median(newdefdat.query('type=="3DM" and deformation=="3D-3D"').cuff_tortuosity)
print(f'Median tortuosity: {tort}')
# plt.title('Histogram of tortuosity values')
# %%
newdef = newdefdat.copy()[newdefdat['deformation'] != "Undeformed"]
newdef['deformation'] = pd.Categorical(newdef['deformation'], categories=['2D-3D', '3D-3D'], ordered=True)
plt.show()
g = sns.catplot(
    kind='violin',
    x='deformation',
    data=newdef.query("fiber_diam in [3,13]"),
    hue='type',
    y='threshold',
    sharey=False,
    errorbar='se',
    palette=pal2d3d,
    col='fiber_diam',
    margin_titles=True,
)
for ax in g.axes.ravel():
    #     ax.set_xlabel('')
    ax.set_ylim([0, None])
g.set_axis_labels('', 'Threshold (mA)')
g.set_titles(col_template="Fiber Diameter: {col_name} μm", row_template="")
# g.fig.set_size_inches(6, 6, forward=True)
# plt.subplots_adjust(top=0.9, wspace=0.4)
# g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# # plt.suptitle('2DEM vs 3DM Thresholds for all samples')
# g.axes.ravel()[0].set_title('3 μm')
# g.axes.ravel()[1].set_title('13 μm')
# plt.show()
g.fig.set_size_inches(12, 5)

# %% dose-response deformed
sns.set(font_scale=2)
sns.set_style('whitegrid')
# could shade in min and max to show response range?
defdr.drop_duplicates().reset_index(drop=True, inplace=True)
plt.figure()
g = sns.relplot(
    kind='line',
    row='fiber_diam',
    data=defdr.query(f"fiber_diam in [3,13] and contact in {main_comparison}"),
    y='percent_activated',
    x='threshold',
    col='deformation',
    hue='nerve_label',
    palette=defpal,
    estimator=None,
    linewidth=3,
    facet_kws={'sharex': False, 'margin_titles': True},
    **{'style': 'modeltype'},
)
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[1:], labels=labs[1:], bbox_to_anchor=[1, 1])

for ax in g.axes.ravel():
    ax.set_ylim([0, 1])
    ax.set_xlim([0, None])
    ax.axhline(0.9, linewidth=2, color='k', linestyle='--', label='saturation')
    ax.axhline(0.1, linewidth=2, color='k', linestyle=':', label='onset')
    # add text at 0.9 and 0.1 for each plot
g.axes[0, 2].text(0.98, 0.84, 'saturation', fontsize=22, transform=g.axes[0][2].transAxes)
g.axes[0, 2].text(0.98, 0.13, 'onset', fontsize=22, transform=g.axes[0][2].transAxes)
g.set_titles(row_template='')

g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion fibers\nactivated (3 μm)')
g.axes[1][0].set_ylabel('Proportion fibers\nactivated (13 μm)')
g._legend.set_title('')
g.set_xlabels('Threshold (mA)')
g.set_titles(col_template="Deformation: {col_name}", row_template='')
g._legend.texts[0].set_text('')
g._legend.texts[5].set_text('')
# %% activation order
newdefdr = defdr.copy()
newdefdr = newdefdr.query("'5R' in nerve_label and deformation!='2D-3D'").sort_values('modeltype')
newdefdr['nerve_label'] = pd.Categorical(newdefdr['nerve_label'], categories=['5R'])
newdefdr['deformation'] = pd.Categorical(newdefdr['deformation'], categories=['Undeformed', '3D-3D'])
sns.set(font_scale=2)
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
# %%
sns.set(font_scale=1.75, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    row='fiber_diam',
    data=newdefdr.query(f"fiber_diam in [3,13] and contact in {main_comparison}"),
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
g.axes[0][0].set_ylabel('Proportion fibers\nactivated (3 μm)')
g.axes[1][0].set_ylabel('Proportion fibers\nactivated (13 μm)')
g.set_xlabels('')
g.legend.remove()
norm = plt.Normalize(0, 1)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=norm)
sm.set_array([])

# Remove the legend and add a colorbar
g.figure.colorbar(
    sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.6, label='Proportion 2DEM fibers activated', pad=0.1
).ax.yaxis.set_ticks_position('left')

# g._legend.set_title('')
# %% recruitment cost:
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
scores = []
datahere = deftomatch.query('deformation in ["3D-3D","Undeformed"]')
for deformation in datahere.deformation.unique():
    for comp in [main_comparison, cath_comparison]:
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
g = sns.FacetGrid(data=scoredat, col='deformation', row='slice', margin_titles=True)
g.map_dataframe(sns.boxplot, y='score2d3d', x='fiber_diam', boxprops={'facecolor': 'white'})
g.map_dataframe(sns.stripplot, y='score2d3d', x='fiber_diam', color='black', size=5, jitter=0.25)
ax.set_xlabel('Fiber Diameter (μm)')
plt.ylabel('AR')
g.set_titles(col_template='{col_name}')
g.set_ylabels('AR')
g.set_xlabels('D (μm)')
# plt.title(f'{deformation} - {comp[0]}')
# %% dose-response onset sat deformed
peses = []
pemeans = []
for stringdat in ['Undeformed', '2D-3D', '3D-3D']:
    subdat = newdefdat.query(f"deformation=='{stringdat}'")
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
    g = sns.FacetGrid(
        compiled_data.query("fiber_diam in [3,13]").rename(
            columns={'fiber_diam': "Fiber Diameter (μm)", 'threshold': "Threshold (mA)"}
        ),
        col="level",
        row="Fiber Diameter (μm)",
        sharey=False,
        margin_titles=True,
    )
    g.map_dataframe(sns.swarmplot, x='type', y='Threshold (mA)', color='black')
    g.map_dataframe(
        sns.lineplot,
        x='type',
        y='Threshold (mA)',
        units='units',
        hue='nerve_label',
        estimator=None,
        color='k',
        palette=defpal,
    )
    plt.subplots_adjust(top=0.9)
    g.set_titles(col_template='{col_name}', row_template='')
    for ax in g.axes.ravel():
        ax.set_xlabel('')
    g.axes[0][0].set_ylabel('Threshold (mA)\n3 μm fibers')
    g.axes[1][0].set_ylabel('Threshold (mA)\n13 μm fibers')
    # plt.gcf().set_size_inches([5, 7])
    plt.legend(title='Sample', bbox_to_anchor=(1.8, 1.5))
    plt.suptitle('Deformation: ' + stringdat)

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
    # plot percent error for each sample DR
    sns.set(font_scale=1.75)
    sns.set_style('whitegrid')
    plt.figure()
    sns.lineplot(data=pes, style='level', y='pe', x='fiber_diam', markers=True, hue='nerve_label')
    plt.xlabel('Fiber Diameter (μm)')
    plt.ylabel('Percent Error')
    # plt.title(f'Percent Error for {stringdat} Onset and Saturation')
    plt.ylim([0, 40])
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

    # plot percent error for population
    sns.set(font_scale=1.75)
    sns.set_style('whitegrid')
    plt.figure()
    sns.lineplot(data=pemean, style='level', y='pe', x='fiber_diam', markers=True)
    plt.legend(title='')
    plt.xlabel('Fiber Diameter (μm)')
    plt.ylabel('Percent Error')
    plt.title(f'Percent Error for {stringdat} pop Onset and Saturation')
    plt.ylim([0, 40])
sns.set(font_scale=2, style='whitegrid')
allpes = pd.concat(peses)
g = sns.relplot(
    data=allpes, kind='line', col='deformation', style='level', y='pe', x='fiber_diam', markers=True, hue='nerve_label'
)
# replace labels
g._legend.texts[0].set_text('')
g._legend.texts[5].set_text('')
g.set_xlabels('Fiber Diameter (μm)')
g.set_titles(col_template='Deformation: {col_name}')
g.set_ylabels("Percent Difference")
plt.xticks([3, 5, 7, 9, 11, 13])
plt.ylim([0, None])
ylimpe = g.axes[0][0].get_ylim()

allpemean = pd.concat(pemeans)
g = sns.relplot(
    data=allpemean, kind='line', col='deformation', style='level', y='pe', x='fiber_diam', markers=True, color='k'
)
g.legend.set_title('')
g.set_xlabels('Fiber Diameter (μm)')
g.set_titles(col_template='Deformation: {col_name}')
g.set_ylabels("Percent Difference")
plt.xticks([3, 5, 7, 9, 11, 13])
plt.ylim(ylimpe)
sns.move_legend(g, 'lower center', bbox_to_anchor=(0.75, 0.5), ncol=1)
# %% deformation mean error threshold
tomatch = newdefdat.query('deformation=="3D-3D"')
defmatched = datamatch(tomatch.query('type=="2DEM"'), tomatch.query('type=="3DM"'), 'threshold').drop(columns='type')
tomatch = newdefdat.query('deformation=="Undeformed"')
ogmatched = datamatch(tomatch.query('type=="2DEM"'), tomatch.query('type=="3DM"'), 'threshold').drop(columns='type')
defsupermatch = pd.concat([defmatched, ogmatched])
# calculate mean, std, and median for each sample
outdata = defmatched.groupby(['nerve_label', 'fiber_diam']).agg(
    {'threshold': ['mean', 'std', 'median'], 'threshold3d': ['mean', 'std', 'median']}
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
plt.ylabel('Percent Difference')
plt.xlabel('Fiber Diam')
plt.title('Mean Threshold Percent Difference\nDeformation=3D-3D')
plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0))
plt.show()
# %% new idea analysis
concats = []
for deftype in ["Undeformed", "3D-3D"]:
    thismatch = deftomatch.query('deformation==@deftype')
    matchednow = datamatch(thismatch.query('type=="2DEM"'), thismatch.query('type=="3DM"'), 'threshold').drop(
        columns='type'
    )
    matchednow['deformation'] = deftype
    concats.append(matchednow)
concats = pd.concat(concats)
concats['deformation'] = pd.Categorical(concats['deformation'], categories=['Undeformed', '3D-3D'], ordered=True)


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


# %%
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

nsimdata = minthresh.query('fiber_diam in [3,13]')  # TODO replace all cath comparison with non
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
g.set_titles('D: {col_name} μm')
g.set_xlabels('2DEM Threshold (mA)')
g.set_ylabels('3DM Threshold (mA)')
# %% LCCC comparison maincath

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
# %% Deformation thresholds unity line
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
# usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
usedata = defsupermatch
usedata['deformation'] = usedata.deformation.astype(str)
nsimdata = usedata.query(
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
    r, p = pearsonr(rdata.threshold, rdata.threshold3d)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    lim = ax.get_xlim()[1]
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} μm: {r ** 2:.2f}')
    ax.set_title(f'{diam} μm')
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='1:1 line')
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])

    # ax.set_yticks(ax.get_xticks())

g.set_xlabels('2DEM Threshold (mA)')
g.set_ylabels('3DM Threshold (mA)')
# %% matching
threeddefmatch = deftomatch.query('deformation=="3D-3D"')
deffinalmatch = datamatch(threeddefmatch.query('type=="2DEM"'), threeddefmatch.query('type=="3DM"'), 'threshold').drop(
    columns='type'
)
# %% Correlation2DEM3DM deformed
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = deffinalmatch.rename(columns={'threshold': 'threshold2DEM'})
for comparison in [
    ['threshold2DEM', 'threshold3d'],
    ['threshold2DEM', 'peri_thk'],
    ['threshold3d', 'peri_thk'],
    # # ['threshold2DEM', 'minimum_efib_distance'],
    ['threshold2DEM', 'peak_second_diff'],
    ['peak_second_diff', 'peri_thk'],
    # # ['peak_second_diff', 'minimum_efib_distance'],
    ['peak_second_z', 'activation_zpos'],
    ['peak_second_diff_node', 'apnode'],
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
# %% BEGIN NEW ANALYSIS
newdata = pd.read_csv("thresh_unmatched_sim3_new.csv")
newdata = addpwfd(newdata, '3')
newdata['contact'] = newdata['sample'].astype(str).str[2].replace({'2': 'cathodic', '0': 'anodic', '1': 'center'})

# %% analyze effect of shifting cuff, datapoints with "sample" ending in 4 or 6
sns.set(style='whitegrid', font_scale=1.75)
shiftdata = newdata[newdata['sample'].astype(str).str[-1].isin(['4', '6'])]
# $rename type from 3D and 2D to 3DM and 2DEM
shiftdata['type'] = shiftdata['type'].replace({'3D': '3DM', '2D': '2DEM'})
ogdata = threshload.query(f'nerve_label=="2L" and contact in {cath_comparison}')
shiftdata = pd.concat([ogdata, shiftdata]).reset_index(drop=True)
# set nerve label as categorical ordered
shiftdata['nerve_label'] = pd.Categorical(shiftdata['nerve_label'], ordered=True, categories=['2Lup', '2L', '2Ldown'])


def fixed_boxplot(*args, label=None, **kwargs):
    sns.boxplot(*args, **kwargs, labels=[label])


# hue_kws = dict(boxprops=dict(edgecolor=color),whiskerprops=dict(color=color),medianprops=dict(color=color),capprops=dict(color=color))
g = sns.catplot(
    data=shiftdata.query('fiber_diam in [3,13]'),
    y='nerve_label',
    x='threshold',
    kind='strip',
    # row='type',
    col='fiber_diam',
    facet_kws=dict(margin_titles=True),
    sharex=False,
    hue='type',
    palette=pal2d3d,
    linewidth=1,
    edgecolor='black',
    jitter=True,
    dodge=True,
    s=50,
)
g.map_dataframe(fixed_boxplot, y='nerve_label', x='threshold', linewidth=3, hue='type', palette=pal2d3d)
# g.map_dataframe(sns.boxplot,data=shiftdata.query(f'type == "{type}"'),y='nerve_label', x='threshold', boxprops=dict(facecolor='None',edgecolor=color),whiskerprops=dict(color=color),medianprops=dict(color=color),capprops=dict(color=color),linewidth=3)
g.set_titles(col_template='D: {col_name} μm', row_template='{row_name}')
g.set_ylabels('')
g.set_xlabels('Threshold (mA)')
for ax in g.axes[:, 0]:
    ax.set_yticks(ax.get_yticks(), ['+4 mm', '+-0mm', '-4 mm'])
# TODO set facecolor as none and edgecolors to the facecolor
# %% compare 2LDSfine to regular 2LDS
sns.set(style='whitegrid', font_scale=1.75)
finedata = newdata[newdata['nerve_label'] == '2Lfine']
unfine_data = threshload[(threshload['nerve_label'] == '2L') & (threshload['contact'] != 'anodic')]
fuf = pd.concat([finedata, unfine_data])
fuf['type'] = fuf['type'].replace({'2D': '2DEM', '3D': '3DM'})
g = sns.catplot(
    data=fuf.query('fiber_diam in [3,13] and type=="3DM"'),
    x='nerve_label',
    y='threshold',
    kind='swarm',
    col='fiber_diam',
    facet_kws=dict(margin_titles=True),
    sharey=False,
    color='w',
    linewidth=1,
    edgecolor='black',
    # jitter=True,
    s=5,
)
g.set_ylabels('Threshold (mA)')
g.set_xlabels('')
g.fig.set_size_inches([12, 5])
# %% imthera notuse
imdata = newdata[newdata['nerve_label'] == '2Lim']
for srcs, name in zip(pd.unique(imdata['active_src_index']), ['monopolar', 'bipolar']):
    thisdata = imdata[imdata['active_src_index'] == srcs]
    g = sns.catplot(
        data=thisdata.query('fiber_diam in [3,13]'),
        y='threshold',
        kind='strip',
        row='fiber_diam',
        facet_kws=dict(margin_titles=True),
        sharey=False,
        hue='type',
        palette=pal2d3d,
        x='type',
        linewidth=1,
        edgecolor='black',
        jitter=True,
        s=50,
    )
    g.set_ylabels('Threshold (mA)')
    g.set_xlabels('')
    g.set_titles(col_template='{col_name}', row_template='{row_name} μm')
    plt.suptitle(name)
# %% compare 2LDSfine to regular 2LDS #TODO figure out why these are different
sns.set(style='whitegrid', font_scale=1.75)
finedata = newdata[newdata['nerve_label'] == '2Lfine']
unfine_data = threshload[(threshload['nerve_label'] == '2L') & (threshload['contact'] != 'anodic')]
fuf = pd.concat([finedata, unfine_data])
fuf['type'] = fuf['type'].replace({'2D': '2DEM', '3D': '3DM'})
g = sns.catplot(
    data=fuf.query('fiber_diam in [3,13] and type=="2DEM" and contact=="cathodic"'),
    x='nerve_label',
    y='threshold',
    kind='strip',
    row='fiber_diam',
    facet_kws=dict(margin_titles=True),
    sharey=False,
    color='w',
    linewidth=1,
    edgecolor='black',
    jitter=True,
    s=50,
)
g.set_ylabels('Threshold (mA)')
g.set_xlabels('')

# %% imthera
imdata = pd.read_csv("thresh_unmatched_sim16_immy.csv")
imdata = addpwfd(imdata, '3')
imdata['contact'] = imdata['sample'].astype(str).str[2].replace({'2': 'cathodic', '0': 'anodic', '1': 'center'})

# inners need to match the center slice
innerid = {}
for inner in pd.unique(imdata.inner):
    innerid[inner] = pd.unique(imdata.query(f'inner == {inner} and type=="2D"')['master_fiber_index'])
    # get all rows where 3D and master fiber index is in the innerid
    imdata.loc[(imdata['type'] == '3D') & (imdata['master_fiber_index'].isin(innerid[inner])), 'inner'] = inner
imdata['type'] = imdata['type'].replace({'2D': '2DEM', '3D': '3DM'})
contact_dict = {'0': 'M1', '1': 'M2', '2': 'C1', '3': 'C2', '4': 'C3', '5': 'C4', '6': 'C5', '7': 'C6'}
imdata['active_src_index'] = imdata['active_src_index'].astype(str).replace(contact_dict)
# %% imthera contact config
sns.set(font_scale=1.75, style='whitegrid')
g = sns.catplot(
    data=imdata.query('fiber_diam in [3,13]'),
    x='active_src_index',
    y='threshold',
    kind='strip',
    row='fiber_diam',
    facet_kws=dict(margin_titles=True),
    sharey='row',
    hue='inner',
    palette='Set1',
    col='type',
    linewidth=0.75,
    s=30,
)
g.set_ylabels('Threshold (mA)')
g.set_xlabels('Contact Configuration')
g.set_titles(col_template='{col_name}', row_template='D: {row_name} μm')
# %%
sns.set(font_scale=1, style='whitegrid')
# contact configurations
ccs = [
    [0, -1, 0, 0, 0, 0],
    [0, 0, 0, 0, -1, 0],
    [-1, 0, 1, 0, 0, 0],
    [0, 0, 0, -1, 0, 1],
    [-1, 1, 0, 0, 0, 0],
    [0, -1, 1, 0, 0, 0],
    [0, 0, 0, -1, 1, 0],
    [0, 0, 0, 0, -1, 1],
]
import matplotlib.pyplot as plt
import numpy as np


def draw_grid(numbers, ax, color=None):
    grid = np.flip(np.array(numbers).reshape(2, 3).T, axis=0)

    cmap = plt.get_cmap("coolwarm", np.max(grid) - np.min(grid) + 2)
    bounds = np.arange(np.min(grid) - 0.5, np.max(grid) + 1.5, 1)
    norm = plt.Normalize(vmin=np.min(grid) - 0.5, vmax=np.max(grid) + 0.5)

    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            # if color:
            #     color=color
            if grid[i, j] == 0:
                color = "black"
            else:
                color = cmap(norm(grid[i, j]))

            circle = plt.Circle((j, i), 0.4, color=color, clip_on=False)
            ax.add_patch(circle)
            ax.text(j, i, grid[i, j], ha='center', va='center', color='white')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim([-0.5, grid.shape[1] - 0.5])
    ax.set_ylim([-0.5, grid.shape[0] - 0.5])
    ax.set_aspect('equal')

    # plt.box(False)


# Example usage
fig, axs = plt.subplots(1, len(ccs))
for i, cc in enumerate(ccs):
    draw_grid(cc, axs[i])
    axs[i].set_title(contact_dict[str(i)])

# make one lone grid with numbers 0,1,2,3,4,5 and no colors
fig, ax = plt.subplots(1, 1)
draw_grid([0, 1, 2, 3, 4, 5], ax, color='black')
ax.set_title('All Contacts')

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
                f'active_src_index == "{contact_config}" and fiber_diam == {fiber_diam} and type == "{t}"'
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
    palette='Set1',
    linewidth=3,
)
plt.ylim(0, 1)
g.set_xlabels('Contact Configuration')
g.set_titles(col_template='D: {col_name} μm')
plt.xticks(range(7))
g.fig.set_size_inches([8, 10])
g.set_titles(row_template='D: {row_name} μm')
# %%
thisdef = deftomatch.query('deformation!="2D-3D" and fiber_diam in [3,13]')
thisdef['deformed'] = thisdef['deformation'] == "3D-3D"

comparecols = ['peak_second_diff', 'peri_thk', 'peri_thk_act_site', 'tortuosity', 'minimum_efib_distance']
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
]
for t, l in zip(leg.get_texts(), new_labs):
    t.set_text(l)
plt.ylim(0, 1)
g.set_ylabels(r'$R^2$', rotation=90)
g.set_titles(col_template="{col_name}", row_template="D: {row_name} μm")

# %% deformation zpos
for stringdat in ['Undeformed', '2D-3D', '3D-3D']:
    sns.set(font_scale=1.75, style='whitegrid')
    newthreshz = newdefdat.query(f'contact in {main_comparison} and deformation=="{stringdat}"')
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
    axs[0].set_title("2DEM-100%")
    axs[1].set_title("3DM-100%")
    plt.suptitle(f'Deformation: {stringdat}', y=1.1)

# %% dose-response
sns.set(font_scale=1.75, style='whitegrid')
defsamples = ["2L", "3R", "5R", "6R"]
und_def_dr = drdat.query(f"nerve_label in {defsamples}")
und_def_dr['nerve_label'] = pd.Categorical(und_def_dr['nerve_label'], categories=defsamples, ordered=True)
plt.figure()
g = sns.relplot(
    kind='line',
    col='nerve_label',
    data=und_def_dr.query(f"fiber_diam in [3] and contact in {main_comparison}"),
    y='percent_activated',
    x='threshold',
    units='nerve_label',
    hue='modeltype',
    palette=pal2d3d,
    estimator=None,
    linewidth=2,
    facet_kws={'sharex': False},
)

g.legend.set_title('')
# sns.move_legend(g,[0.5,0.5])
for ax in g.axes.ravel():
    ax.set_xlabel('Threshold (mA)')
    ax.set_ylim([0, 1])
    ax.set_xlim([0, None])
g.set_titles(row_template='', col_template='{col_name}')

# g.axes[0][0].set_xlabel('')
g.set_ylabels('Proportion activated')
# g.axes[0][0].axhline(0.9, linewidth=2, color='k', linestyle='-', label='saturation')
# g.axes[0][0].axhline(0.1, linewidth=2, color='k', linestyle='--', label='onset')
# g.axes[1][0].axhline(0.9, linewidth=2, color='k', linestyle='-', label='saturation')
# g.axes[1][0].axhline(0.1, linewidth=2, color='k', linestyle='--', label='onset')
# # add text at 0.9 and 0.1 for each plot
# g.axes[0][0].text(1.0, 0.87, 'saturation', fontsize=18, transform=g.axes[0][0].transAxes)
# g.axes[0][0].text(1.0, 0.1, 'onset', fontsize=18, transform=g.axes[0][0].transAxes)
# g.axes[1][0].text(1.0, 0.87, 'saturation', fontsize=18, transform=g.axes[1][0].transAxes)
# g.axes[1][0].text(1.0, 0.1, 'onset', fontsize=18, transform=g.axes[1][0].transAxes)
# g.fig.set_size_inches([7, 10])

# %% dose-response
sns.set(font_scale=1.75, style='whitegrid')
defsamples = ["2L", "3R", "5R", "6R"]
newdef = defdr.copy()[defdr['deformation'] != "Undeformed"]
newdef['deformation'] = pd.Categorical(newdef['deformation'], categories=['2D-3D', '3D-3D'], ordered=True)
plt.figure()
g = sns.relplot(
    kind='line',
    col='nerve_label',
    style='deformation',
    data=newdef.query(
        f"fiber_diam in [3] and contact in {main_comparison} and nerve_label in {defsamples} and deformation!='Undeformed'"
    ),
    y='percent_activated',
    x='threshold',
    units='nerve_label',
    hue='modeltype',
    palette=pal2d3d,
    estimator=None,
    linewidth=2,
    facet_kws={'sharex': False, 'margin_titles': True},
)

g.legend.set_title('')
# sns.move_legend(g,[0.5,0.5])
for ax in g.axes.ravel():
    ax.set_ylim([0, 1])
    ax.set_xlim([0, None])
g.set_titles(row_template='', col_template='{col_name}')

# g.axes[0][0].set_xlabel('')
g.set_ylabels('Proportion activated')
g.set_xlabels('Threshold (mA)')
# g.axes[0][0].axhline(0.9, linewidth=2, color='k', linestyle='-', label='saturation')
# g.axes[0][0].axhline(0.1, linewidth=2, color='k', linestyle='--', label='onset')
# g.axes[1][0].axhline(0.9, linewidth=2, color='k', linestyle='-', label='saturation')
# g.axes[1][0].axhline(0.1, linewidth=2, color='k', linestyle='--', label='onset')
# # add text at 0.9 and 0.1 for each plot
# g.axes[0][0].text(1.0, 0.87, 'saturation', fontsize=18, transform=g.axes[0][0].transAxes)
# g.axes[0][0].text(1.0, 0.1, 'onset', fontsize=18, transform=g.axes[0][0].transAxes)
# g.axes[1][0].text(1.0, 0.87, 'saturation', fontsize=18, transform=g.axes[1][0].transAxes)
# g.axes[1][0].text(1.0, 0.1, 'onset', fontsize=18, transform=g.axes[1][0].transAxes)
# g.fig.set_size_inches([7, 10])
# %% dose-response deformed vs undeformed 3D only
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
# could shade in min and max to show response range?
defdefdr = defdr.query("modeltype=='3DM' and deformation!='3D-3D'")
defdefdr['deformed'] = defdefdr['deformation'] == "2D-3D"

plt.figure()
g = sns.relplot(
    kind='line',
    data=defdefdr.query("fiber_diam in [3]"),
    y='percent_activated',
    x='threshold',
    # hue='nerve_label',
    col="nerve_label",
    # palette=defpal,
    estimator=None,
    linewidth=3,
    color=pal2d3d[1],
    facet_kws={'sharex': False, 'margin_titles': True},
    **{'style': 'deformed'},
)

# sns.move_legend(g, 'lower center', bbox_to_anchor=(0.9, 0.4), ncol=1)
g.axes[0][0].set_xlabel('')
g.set_ylabels('Proportion Activated')
g.set_xlabels('Threshold (mA)')
g.set_titles(col_template="{col_name}", row_template='{row_name}')
