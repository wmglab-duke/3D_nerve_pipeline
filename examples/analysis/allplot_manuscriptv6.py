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
        data.loc[data['nsim'] == nsim, 'pulse_width'] = pulse_width
        data.loc[data['nsim'] == nsim, 'fiber_diam'] = fiber_diam
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

# gogo = "deformedasc"

#%% set which comparison to run
# base data
threshload = pd.read_csv('thresh_unmatched_sim3.csv').query('sim==3')
# set all inners and outers where 'type' is 3D to 0
threshload.loc[(threshload['type'] == '3D'), 'inner'] = 0
threshload.loc[(threshload['type'] == '3D'), 'outer'] = 0

# inners need to match the cathodic leading contact
for i, row in threshload.iterrows():
    if row['type'] == '3D':
        # quickinner is only where the sample ends in a "2", must convert int to string
        quickinner = threshload[threshload['sample'].astype(str).str[-1] == '2']
        nerver = row.nerve_label[:2]
        data = quickinner.query(
            f"nerve_label == '{nerver}' and master_fiber_index == @row.master_fiber_index and type == '2D' and nsim==0"
        )
        assert len(data) == 1
        threshload.loc[i, 'inner'] = data['inner'].values[0]
        threshload.loc[i, 'outer'] = data['outer'].values[0]
        threshload.loc[i, 'fiber'] = data['fiber'].values[0]

# threshload.to_csv('thresh_unmatched_sim3.csv')
if gogo == "initial":  # remove all where nerve_label length is >2
    threshload = threshload[threshload['nerve_label'].str.len() < 3]
# elif gogo=="deformedasc": #remove all where nerve_label does not contain "asc"
#     threshload = threshload[threshload['nerve_label'].str.contains("asc")]
#     #check the third digit of the sample number is 2, in that case, "contact" is cathodic. If 0, anodic
#     threshload['contact'] = threshload['sample'].astype(str).str[2].replace({'2': 'cathodic', '0': 'anodic'})
#%%Setup
threshload['type'] = threshload['type'].replace({'2D': '2DEM', '3D': '3DM'})
threshload = addpwfd(threshload, '3')
threshload['fiber_diam'] = threshload['fiber_diam'].astype(int)
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

pal2d3d = ['#d95f02', '#7570b3']

# sys.exit('prepdone')
#%% sample thresholds with unity line
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
# usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
usedata = matched
nsimdata = usedata.query('fiber_diam in [3,13]')
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
    print(f'{diam} μm: {r**2:.2f}')
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
#%% threshold violinplot
# generate boxplot of 2DEM and 3DM thresholds
# sns.set(font_scale=1.75)
sns.set_style('whitegrid')
g = sns.catplot(
    kind='violin',
    col='fiber_diam',
    row='nerve_label',
    data=threshload,
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
# g.fig.set_size_inches(5, 5, forward=True)
plt.subplots_adjust(top=0.9, wspace=0.4)
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# plt.suptitle('2DEM vs 3DM Thresholds for all samples')
g.axes.ravel()[0].set_title('3 μm')
g.axes.ravel()[1].set_title('13 μm')
plt.show()
# calculate mean, std, and median for each sample
outdata = matched.groupby(['nerve_label', 'fiber_diam']).agg(
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
plt.title('Mean Threshold Percent Difference')
plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0))
plt.show()
sns.lineplot(data=outdata, x='fiber_diam', y='std_diff', hue='nerve_label', palette='colorblind', marker='o')
plt.ylabel('Percent Difference')
plt.xlabel('Fiber Diam')
plt.title('Standard Deviation Percent Difference')
plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0))
plt.show()
sns.lineplot(data=outdata, x='fiber_diam', y='median_diff', hue='nerve_label', palette='colorblind', marker='o')
plt.ylabel('Percent Difference')
plt.xlabel('Fiber Diam')
plt.title('Median Threshold Percent Difference')
plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0))
plt.show()
#%% threshold vals
for nsim in [0, 5]:
    dadata = matched.query(f'nsim=={nsim}')
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
#%% threshold variances and coefficient of variation intrafascicle and inter
sns.set(font_scale=1.75, style='whitegrid')
vardat = repeated.query('sim==3 and contact=="cathodic"')
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
#%% dose-response
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
                    f'nerve_label=="{sampname}" and fiber_diam=={fiber_diam} and contact!="anodic" and modeltype=="{modeltype}"'
                )
                thisdat = thisdat.sort_values(by='threshold')
                thisdat = thisdat.reset_index(drop=True)
                thisdat = thisdat.loc[thisdat.percent_activated >= percent / 100]
                thesevals.append(thisdat.threshold.iloc[0])
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
    data=drdat.query("fiber_diam in [3,13] and contact != 'anodic'"),
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
#%% all dose-response
newdr = drdat.copy()
sns.set(font_scale=2)
sns.set_style('whitegrid')
newdr['percent_activated2d'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for i, row in newdr.iterrows():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdr.query(
        'modeltype == "2DEM" and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim and master_fiber_index == @row.master_fiber_index and contact == "cathodic"'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    newdr.loc[i, 'percent_activated2d'] = val
#%%
# first find the mean threshold across all inners
grouped = newdr.query('modeltype == "2DEM" and contact == "cathodic"').groupby(
    ['nerve_label', 'fiber_diam', 'sim', 'inner']
)
meanthresh = grouped.threshold.mean().reset_index().dropna()
# add to meanthresh the percentile of each in`ner's threshold
meanthresh['percentile'] = np.nan
for i, row in meanthresh.iterrows():
    thisdat = meanthresh.query('nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim')
    # percentile is number of mean thresholds <= this inner's threshold / total number of inners
    meanthresh.loc[i, 'percentile'] = len(thisdat.query('threshold <= @row.threshold')) / len(thisdat)
newdr['inner_activation_2D_percent'] = np.nan
analyze = newdr.copy().query("contact!='anodic'")
for i, row in analyze.iterrows():
    thisdat = meanthresh.query('nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim')
    # find the inner percentile for this fiber
    val = thisdat.query('inner == @row.inner').percentile.values[0]
    assert not val is np.nan
    analyze.loc[i, 'inner_activation_2D_percent'] = val

# could shade in min and max to show response range?
analyze['nerve_label'] = pd.Categorical(
    analyze['nerve_label'], ordered=True, categories=['2L', '2R', '3R', '5R', '6L', '6R']
)
#%%
sns.set(font_scale=2, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    row='fiber_diam',
    data=newdr.query("fiber_diam in [3,13] and contact != 'anodic'"),
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
)

g.set_titles(row_template='', col_template='{col_name}')
g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion of fibers activated\nFiber diameter: 3 μm')
g.axes[1][0].set_ylabel('Proportion of fibers activated\nFiber diameter: 13 μm')
g.set_xlabels('')
g._legend.set_title('Percentage of\n2D fibers activated')
#%% organization compare type with recruitment order
# remove all samples which dont end with a 2
sns.reset_orig()
sns.set(font_scale=1.25, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
# threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]
# # subtract 1 from sample if type is 3DM
# threshdat['sample'] = threshdat.apply(lambda x: x['sample'] - 1 if x.type == "3DM" else x['sample'], axis=1)
data = newdr.query("fiber_diam in [3,13] and contact != 'anodic' and nerve_label=='5R'")
# apply minmax normalization within sample and nsim
g = sns.FacetGrid(
    data,
    col="nerve_label",
    row='fiber_diam',
    sharey=False,
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x='modeltype',
    y='percent_activated',
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
    y='percent_activated',
    palette='rainbow',
    hue='inner',
)

plt.subplots_adjust(top=0.9)
g.fig.set_size_inches(4, 6)
g.set_titles(col_template='{col_name}', row_template='')
g.set_xlabels('')
g.axes[0, 0].set_ylabel('Percent\nFiber Diam: 3μm')
g.axes[1, 0].set_ylabel('Percent\nFiber Diam: 13μm')

# handles, labs = plt.gca().get_legend_handles_labels()
# plt.legend(title='inner', handles=handles[14:], labels=labs[14:], bbox_to_anchor=[1, 2])


# plt.suptitle('thresholds compared between 2DEM (cathodic) and 3DM')
g = sns.FacetGrid(
    data,
    col="nerve_label",
    row='fiber_diam',
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
g.fig.set_size_inches(3, 6)
g.set_titles(col_template='{col_name}', row_template='')
g.set_xlabels('')

# handles, labs = plt.gca().get_legend_handles_labels()
# plt.legend(title='inner', handles=handles[14:], labels=labs[14:], bbox_to_anchor=[1, 2])

g.axes[0, 0].set_ylabel('Threshold (mA)\nD: 3μm')
g.axes[1, 0].set_ylabel('Threshold (mA)\nD: 13μm')
# plt.suptitle('thresholds compared between 2DEM (cathodic) and 3DM')
#%% 3d activation as a percent of 2d activation
sns.set(font_scale=2)
sns.set_style('whitegrid')
# go through drdat and for each threshold, find the percentage of fibers of the opposite modeltype that are activated
drdat['opposite_model_activation'] = np.nan
for i, row in drdat.iterrows():
    # need to remove all anodic data
    thisdat = drdat.query(
        'modeltype != @row.modeltype and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim'
    )
    # percent is number of thresholds less than or equal to this threshold divided by total number of thresholds
    drdat.loc[i, 'opposite_model_activation'] = len(thisdat.query('threshold <= @row.threshold')) / len(thisdat)
#%%final zpos
sns.set(font_scale=1.75, style='whitegrid')
newthreshz = threshload.copy()
newthreshz['activation_zpos'] = newthreshz['activation_zpos'] / 10000
# plt.suptitle('Activation z-position for 3um fibers')
g = sns.displot(
    data=newthreshz.query("nsim in [0,5]").rename(columns={'nerve_label': 'Sample'}),
    y='activation_zpos',
    col='type',
    row='fiber_diam',
    hue='Sample',
    facet_kws={'sharex': False},
    kind='kde',
    palette='colorblind',
    common_norm=False,
    legend=False,
)
ylim = g.axes[0][0].get_ylim()
for ax in g.axes.ravel():
    ax.set_title('')
for i, diam in enumerate([3, 13]):
    g.axes[i][0].set_xlim(reversed(g.axes[i][0].get_xlim()))
    g.axes[i][0].set_ylabel(f'Activation Location (cm)\n{diam} μm fibers')
g.fig.set_size_inches([6, 8])
zloadmono = pd.read_csv('thresh_unmatched_sim10.csv')
zloadmono = addpwfd(zloadmono, '10')
zloadmono['fiber_diam'] = zloadmono['fiber_diam'].astype(int)
zloadmono['activation_zpos'] = zloadmono['activation_zpos'] / 10000
g = sns.displot(
    data=zloadmono.query("nsim in [0,5]").rename(columns={'nerve_label': 'Sample'}),
    y='activation_zpos',
    col='type',
    row='fiber_diam',
    hue='Sample',
    facet_kws={'sharex': False},
    kind='kde',
    palette='colorblind',
    common_norm=False,
    legend=False,
)
for ax in g.axes.ravel():
    plt.ylim(ylim)
    ax.set_title('')
for i, diam in enumerate([3, 13]):
    g.axes[i][0].set_xlim(reversed(g.axes[i][0].get_xlim()))
    g.axes[i][0].set_ylabel(f'Activation Location (cm)\n{diam} μm fibers')
g.fig.set_size_inches([6, 8])

#%% Correlation2DEM3DM
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = matched.rename(columns={'threshold': 'threshold2DEM'})
for comparison in [
    ['threshold2DEM', 'threshold3d'],
    ['threshold2DEM', 'peri_thk'],
    ['threshold3d', 'peri_thk'],
]:
    corrs = usedata.groupby(['sample', 'fiber_diam', 'contact', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
    corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    corrs['contact'] = pd.Categorical(corrs['contact'], categories=['cathodic', 'anodic'], ordered=True)
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
#%%correlations monopolar new
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
usedata['type'] = usedata['type'].replace({'2D': '2DEM', '3D': '3DM'})
usedata = datamatch(usedata.query('type=="2DEM"'), usedata.query('type=="3DM"'), 'threshold').drop(columns='type')
# add new EMsample column to matched dataframe, which is the nerve label plus the first letter of the contact type capitalized
usedata['EMsample'] = usedata['nerve_label'] + usedata['contact'].str[0].str.upper()
usedata = usedata.rename(columns={'threshold': 'threshold2DEM'})
for comparison in [
    ['threshold2DEM', 'threshold3d'],
    ['threshold2DEM', 'peri_thk'],
    ['threshold3d', 'peri_thk'],
    # ['threshold3d', 'peri_thk_act_site'],
]:  # TODO: make this into a function
    corrs = usedata.groupby(['sample', 'fiber_diam', 'contact', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
    corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    corrs['contact'] = pd.Categorical(corrs['contact'], categories=['cathodic', 'anodic'], ordered=True)
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
#%%
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
#%% tortuosity lineplot
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
nsimdata = threshload.query('type=="3DM"')
corrs = (
    nsimdata.groupby(['sample', 'fiber_diam', 'nerve_label', 'type'])['threshold', 'tortuosity'].corr().iloc[0::2, -1]
)
corrs = corrs.reset_index().rename(columns={'threshold': 'correlation'})
corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
# print(corrs)
for fiber_diam in pd.unique(corrs['fiber_diam']):
    for nerve in pd.unique(nsimdata['nerve_label']):
        # calculate slope of linear regression of threshold vs tortuosity
        this_data = nsimdata.query(f'fiber_diam=={fiber_diam} and nerve_label=="{nerve}"')
        this_data['threshold'] = this_data['threshold'] / np.amax(this_data['threshold'])
        slope, intercept, r_value, p_value, std_err = stats.linregress(this_data[['threshold', 'tortuosity']])
        # add slope to corrs where fiber_diam and nerve_label match
        corrs.loc[(corrs['fiber_diam'] == fiber_diam) & (corrs['nerve_label'] == nerve), 'slope'] = slope
plt.gcf().set_size_inches(4, 4)
sns.swarmplot(data=corrs, y='slope', x='nerve_label', hue='fiber_diam', palette='RdPu')
plt.legend(bbox_to_anchor=[1, 1], title='D (μm)')
# plt.yscale('log')
plt.ylabel(r'slope $(\frac{mA/mA}{μm/μm})$')
plt.xlabel('')
corrs['R2'] = corrs['tortuosity'] ** 2
means = corrs.groupby(['type', 'fiber_diam']).agg(np.mean)
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
for ax, r in zip(g.axes.ravel(), means['tortuosity']):
    # ax.text(0.5, 0.9, f'\nmean r={r:.2f}', transform=ax.transAxes)
    print(r"$R^2$ = " + f"{r:.2f}")
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# plt.xlabel('3DM Threshold (mA)')
g.set_titles(col_template='{col_name} μm')
g.fig.set_size_inches([10, 5])
plt.show()
sys.exit()
#%% Percent Error
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

#%% Dose-response info all diams
plt.figure()
levels = {
    'onset': 10,
    'saturation': 90,
}
grouped = threshload.query('sim==3').groupby(['sample', 'fiber_diam', 'sim', 'type', 'nerve_label', 'model', 'nsim'])
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


#%% organization compare type
# remove all samples which dont end with a 2
sns.reset_orig()
sns.set(font_scale=1.75, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]
# subtract 1 from sample if type is 3DM
threshdat['sample'] = threshdat.apply(lambda x: x['sample'] - 1 if x.type == "3DM" else x['sample'], axis=1)
data = threshdat.query('nsim in [0,5] and sim == 3')
# apply minmax normalization within sample and nsim
g = sns.FacetGrid(
    data.rename(columns={"threshold": "Threshold (mA)", "nerve_label": "Sample"}),
    col="Sample",
    row='fiber_diam',
    sharey=False,
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x='type',
    y='Threshold (mA)',
    units='master_fiber_index',
    estimator=None,
    color='k',
    linewidth=1,
    alpha=0.25,
)
g.map_dataframe(sns.stripplot, jitter=False, linewidth=1, x='type', y='Threshold (mA)', palette='binary')

plt.subplots_adjust(top=0.9)
g.fig.set_size_inches(10, 6)
g.set_titles(col_template='{col_name}', row_template='')
g.axes.ravel()[0].set_ylabel('Threshold (mA)\nFiber Diam: 3μm')
g.axes.ravel()[6].set_ylabel('Threshold (mA)\nFiber Diam: 13μm')
# plt.suptitle('thresholds compared between 2DEM (cathodic) and 3DM')
#%% organization compare nsim
# remove all samples which dont end with a 2
sns.reset_orig()
sns.set(font_scale=1.75, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]
data = threshdat.query('nsim in [0,5] and sim == 3')
# apply minmax normalization within sample and nsim
data['Threshold (mA)'] = data.groupby(['sample', 'nsim', 'type'])['threshold'].transform(
    lambda x: (x - x.min()) / (x.max() - x.min())
)
# set nsim as categorical
data['Fiber Diameter'] = data['nsim'].astype('category')
# nsim 0 is three and nsim 5 is thirteen
data['Fiber Diameter'].cat.rename_categories({0: '3', 5: '13'}, inplace=True)
g = sns.FacetGrid(data, col="nerve_label", row='type', sharey=False, margin_titles=True)
g.map_dataframe(
    sns.lineplot,
    x='Fiber Diameter',
    y='Threshold (mA)',
    units='master_fiber_index',
    estimator=None,
    color='k',
    alpha=0.25,
)
g.map_dataframe(sns.stripplot, linewidth=1, x='Fiber Diameter', y='Threshold (mA)', palette='binary', jitter=False)
# plt.subplots_adjust(top=0.9)
# plt.suptitle('min-max normalized thresholds compared between 3 um and 13 um thresholds (cath 2DEMEM)')
g.set_titles(col_template='{col_name}', row_template='')
plt.savefig('matchsim.png', dpi=400)
g.axes.ravel()[0].set_ylabel('Threshold (mA)\nModel: 2DEM')
g.axes.ravel()[6].set_ylabel('Threshold (mA)\nModel: 3DM')
plt.gcf().set_size_inches([10, 6])
#%% stimfigs
fig, axs = plt.subplots(2, 2)
axs[0][0].plot([0, 1, 1, 2, 2, 3, 3, 4], [0, 0, -1, -1, 1, 1, 0, 0], 'k', linewidth=3)
axs[1][0].plot([0, 1, 1, 2, 2, 3, 3, 4], [0, 0, 1, 1, -1, -1, 0, 0], 'k', linewidth=3)
axs[0][1].plot([0, 1, 1, 2, 2, 3, 3, 4], [0, 0, -1, -1, 0, 0, 0, 0], 'k', linewidth=3)
axs[0][1].set_ylim(axs[0][0].get_ylim())
axs[1][1].plot([0, 1, 1, 2, 2, 3, 3, 4], [0, 0, 0, 0, 0, 0, 0, 0], 'k', linewidth=3)
axs[1][1].set_ylim(axs[0][0].get_ylim())
for ax in axs.ravel():
    ax.axis('off')
#%% organizatino compare nsim
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]

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
ax = sns.boxplot(data=scoredat, x='score2d3d', y='fiber_diam', boxprops={'facecolor': 'white'})
ax.set_ylabel('Fiber Diameter (μm)')
plt.xlabel('RC')
plt.gcf().set_size_inches([3, 5])
#%% organizatino compare type
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]

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
ax = sns.boxplot(data=scoredat, x='scoresmolbeeg', y='type', boxprops={'facecolor': 'white'})
plt.xlabel('RC')
plt.ylabel('')
plt.gcf().set_size_inches([3, 5])
#%% peri-thk with 2DEM3DM
sns.set(font_scale=2)
sns.set_style('whitegrid')
nsimdata = threshload.query('fiber_diam in [3,13]')
actlocdata = nsimdata.query('type=="3DM"')
actlocdata['type'] = '3DM*'
actlocdata['peri_thk'] = actlocdata['peri_thk_act_site']
nsimdata = pd.concat([nsimdata, actlocdata])
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
g.axes.flat[3].set_xlabel('Perineurium Thickness (μm)')
g.axes.flat[4].set_xlabel('Perineurium Thickness (μm)')
g.axes.flat[5].set_xlabel('Perineurium Thickness (μm)')
g.axes.flat[0].set_ylabel('Threshold (mA)')
g.axes.flat[3].set_ylabel('Threshold (mA)')

for ax, r, letter in zip(
    g.axes.ravel(order='F'), means['peri_thk'], ['A', '', 'B', '', 'C', '']
):  # TODO: these are not ordering correctly
    nl = '\n'
    ax.text(0, 0.9, fr'{letter}{nl}$\bar{{r}}$={r:.2f}', transform=ax.transAxes)
# plt.xlabel('3DM Threshold (mA)')
g.set_titles(col_template='{col_name}', row_template='Fiber Diameter = {row_name}')
# g.fig.set_size_inches(10, 6)
plt.show()
#%%newcorr
supercorr = []
sns.reset_orig()
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = matched.rename(columns={'threshold': 'threshold2DEM'})
comparison = ['threshold2DEM', 'threshold3d']
corrs = usedata.groupby(['sample', 'fiber_diam', 'contact', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
corrs['contact'] = pd.Categorical(corrs['contact'], categories=['cathodic', 'anodic'], ordered=True)
plt.figure()
corrs.rename(inplace=True, columns={"nerve_label": "Sample"})
sns.scatterplot(data=corrs, x='fiber_diam', y='correlation', hue='Sample', s=100, palette='colorblind', style='contact')
ax = sns.lineplot(
    data=corrs,
    x='fiber_diam',
    y='correlation',
    style='contact',
    hue='Sample',
    legend=False,
    palette='colorblind',
)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.gca().set_ylim([-1, 1])
plt.xlabel('Fiber Diameter (μm)')
plt.title(f'Correlation between {comparison[0]} and {comparison[1]}', pad=20)
plt.figure()
# g = sns.FacetGrid(data=corrs,col='contact')
# g.map_dataframe(sns.stripplot, hue='Sample', y='correlation',dodge=True, palette='colorblind')
# g.map_dataframe(sns.boxplot,y='correlation', boxprops={'facecolor':'None'},whis=100)
sns.stripplot(data=corrs, x='contact', hue='Sample', y='correlation', dodge=True, palette='colorblind')
sns.boxplot(data=corrs, x='contact', y='correlation', boxprops={'facecolor': 'None'}, whis=100)
# plt.subplots_adjust(top=0.8)
plt.title(f'Correlation between \n{comparison[0]} and {comparison[1]}', pad=20)
plt.gca().set_ylim([-1, 1])
sns.move_legend(plt.gca(), "upper left", bbox_to_anchor=(1, 1))
corrs['pole'] = 'bi'
corrs['nerve_label'] = corrs['Sample']
supercorr.append(corrs)

sns.set(font_scale=1.75)
sns.set_style('whitegrid')
usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
usedata['type'] = usedata['type'].replace({'2D': '2DEM', '3D': '3DM'})
usedata = datamatch(usedata.query('type=="2DEM"'), usedata.query('type=="3DM"'), 'threshold').drop(columns='type')
# add new EMsample column to matched dataframe, which is the nerve label plus the first letter of the contact type capitalized
usedata['EMsample'] = usedata['nerve_label'] + usedata['contact'].str[0].str.upper()
usedata = usedata.rename(columns={'threshold': 'threshold2DEM'})
comparison = ['threshold2DEM', 'threshold3d']  # TODO: make this into a function
corrs = usedata.groupby(['sample', 'fiber_diam', 'contact', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
corrs['contact'] = pd.Categorical(corrs['contact'], categories=['cathodic', 'anodic'], ordered=True)
plt.figure()
sns.scatterplot(
    data=corrs, x='fiber_diam', y='correlation', hue='nerve_label', s=100, palette='colorblind', style='contact'
)
ax = sns.lineplot(
    data=corrs,
    x='fiber_diam',
    y='correlation',
    style='contact',
    hue='nerve_label',
    legend=False,
    palette='colorblind',
)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.gca().set_ylim([-1, 1])
plt.title(f'Correlation between {comparison[0]} and {comparison[1]} (MM)', pad=20)
plt.figure()
# g = sns.FacetGrid(data=corrs,col='contact')
# g.map_dataframe(sns.stripplot, hue='Sample', y='correlation',dodge=True, palette='colorblind')
# g.map_dataframe(sns.boxplot,y='correlation', boxprops={'facecolor':'None'},whis=100)
sns.stripplot(data=corrs, x='contact', hue='nerve_label', y='correlation', dodge=True, palette='colorblind')
sns.boxplot(data=corrs, x='contact', y='correlation', boxprops={'facecolor': 'None'}, whis=100)
# plt.subplots_adjust(top=0.8)
plt.title(f'Correlation between \n{comparison[0]} and {comparison[1]} (MM)', pad=20)
plt.gca().set_ylim([-1, 1])
sns.move_legend(plt.gca(), "upper left", bbox_to_anchor=(1, 1))
corrs['pole'] = 'mono'
supercorr.append(corrs)

sns.catplot(
    kind='strip', col='contact', data=pd.concat(supercorr), hue="pole", x="nerve_label", y="correlation", dodge=True
)
#%% BEGIN DEFORMATION ANALYSIS
# match fascicle thresholds
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshes = pd.read_csv('thresh_unmatched_sim3.csv').query('sim==3')
# remove all rows where nerve label contains "asc" and "sample" contains "3"
threshes = threshes[~(threshes['nerve_label'].str.contains('asc') & threshes['sample'].astype(str).str.contains('3'))]
# duplicate all samples where nerve label contains "def" and "sample" contains "3"
defdat = threshes[threshes['nerve_label'].str.contains('def') & threshes['sample'].astype(str).str.contains('3')]
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
newdefdat['contact'] = newdefdat['sample'].astype(str).str[2].replace({'2': 'cathodic', '0': 'anodic'})
# remove all aodic
newdefdat = newdefdat[newdefdat['contact'] != 'anodic']
# set deformation as ordered categorical
newdefdat['deformation'] = pd.Categorical(
    newdefdat['deformation'], categories=['Undeformed', '2D-3D', '3D-3D'], ordered=True
)
newdefdat['nerve_label'] = pd.Categorical(newdefdat['nerve_label'], categories=['2L', '3R', '5R', '6R'], ordered=True)
# remove unused colors from palette
defpal = [sns.color_palette('colorblind')[ind] for ind in [0, 2, 3, 5]]
defdefcomp = newdefdat.query('type=="3DM" and deformation != "2D-3D"')
defdefcomp['deformed'] = defdefcomp['deformation'] != "Undeformed"
#%%
sns.set(font_scale=2, style='whitegrid')
g = sns.catplot(
    kind='violin',
    col='deformation',
    data=newdefdat.query("nsim in [0,5]"),
    x='nerve_label',
    y='threshold',
    sharey=False,
    errorbar='se',
    palette=pal2d3d,
    row='fiber_diam',
    hue='type',
    margin_titles=True,
)
for ax in g.axes.ravel():
    # ax.set_xlabel('')
    ax.set_ylim([0, None])
g.set_axis_labels('', 'Threshold (mA)')
g.set_titles(col_template="Deformation: {col_name}", row_template="Fiber Diameter: {row_name} μm")
g.legend.set_title('')

# g.fig.set_size_inches(6, 6, forward=True)
# plt.subplots_adjust(top=0.9, wspace=0.4)
# g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# # plt.suptitle('2DEM vs 3DM Thresholds for all samples')
# g.axes.ravel()[0].set_title('3 μm')
# g.axes.ravel()[1].set_title('13 μm')
plt.show()
g = sns.catplot(
    kind='violin',
    col='deformation',
    data=newdefdat.query("nsim in [0,5]"),
    x='type',
    y='threshold',
    sharey=False,
    errorbar='se',
    palette=pal2d3d,
    row='fiber_diam',
    margin_titles=True,
)
for ax in g.axes.ravel():
    #     ax.set_xlabel('')
    ax.set_ylim([0, None])
g.set_axis_labels('', 'Threshold (mA)')
g.set_titles(col_template="Deformation: {col_name}", row_template="Fiber Diameter: {row_name} μm")
# g.fig.set_size_inches(6, 6, forward=True)
# plt.subplots_adjust(top=0.9, wspace=0.4)
# g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# # plt.suptitle('2DEM vs 3DM Thresholds for all samples')
# g.axes.ravel()[0].set_title('3 μm')
# g.axes.ravel()[1].set_title('13 μm')
# plt.show()
#%%fascicle analysis
# matchednon = datamatch(newdefdat.query('type=="2DEM" and deformed==False'), newdefdat.query('type=="3DM" and deformed==False'), 'threshold').drop(columns='type')
# matcheddef = datamatch(newdefdat.query('type=="2DEM" and deformed==True'), newdefdat.query('type=="3DM" and deformed==True'), 'threshold').drop(columns='type')
# #calculate a new column which is the mean threshold for that inner
# allmatch = pd.concat([matchednon, matcheddef])
# allmatch['mean_thresh'] = allmatch.groupby(['sample', 'nerve_label', 'inner'])['threshold'].transform('mean')
# first, get only the 3D data
newmatch = newdefdat.query('type=="3DM"')
# for all 2DEM data, aggregate the mean threshold across inners
new2dem = (
    newdefdat.query('type=="2DEM"')
    .groupby(['sample', 'nerve_label', 'fiber_diam', 'inner', 'deformation'])['threshold']
    .mean()
    .reset_index()
)
# drop all rows where threshold is nan
new2dem = new2dem.dropna()
# rename newmatch threshold to threshold3d
newmatch = newmatch.rename(columns={'threshold': 'threshold3d'})
# merge newmatch and new2dem
allmatch = pd.merge(newmatch, new2dem, on=['nerve_label', 'fiber_diam', 'deformation', 'inner'])
allmatch = allmatch.rename(columns={'threshold': 'mean_thresh'})
#%%
sns.set(font_scale=2)
sns.set_style('whitegrid')
# usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
usedata = allmatch.copy()
nsimdata = usedata.query('fiber_diam in [3,13]')
g = sns.relplot(
    data=nsimdata.rename(columns={'nerve_label': 'Sample'}),
    kind='scatter',
    col='deformation',
    x='mean_thresh',
    y='threshold3d',
    hue='Sample',
    s=60,
    palette=defpal,
    facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
    row='fiber_diam',
)

for diam, defstatus, ax in zip(
    [3, 3, 3, 13, 13, 13], ('Undeformed', '2D-3D', '3D-3D', 'Undeformed', '2D-3D', '3D-3D'), g.axes.ravel()
):
    # ax.set_title(f'{diam}, {defstatus}')
    limmax = max([ax.get_xlim()[1], ax.get_ylim()[1]])
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation == "{defstatus}"')
    r, p = pearsonr(rdata.mean_thresh, rdata.threshold3d)
    perc = sum(rdata.mean_thresh > rdata.threshold3d) / len(rdata.mean_thresh)
    # add correlation to plot
    ax.text(0.02, 0.9, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    ax.set_title(f'Fiber Diameter: {diam} μm')
    ax.plot([0, rdata.mean_thresh.max()], [0, rdata.mean_thresh.max()], '--k', linewidth=2, label='1:1 line')
    # ax.set_xlabel('2DEM Threshold (mA)')
    # ax.set_yticks(ax.get_xticks())
    ax.set_xlim([0, limmax])
    ax.set_ylim([0, limmax])
plt.legend()
g.set_axis_labels('Mean 2DEM Threshold\n(mA, by fascicle)', '3DM Threshold\n(mA, by fiber)')
g.set_titles(col_template="Deformation: {col_name}", row_template="Fiber Diameter: {row_name} μm")

# g.axes.ravel()[0].set_ylabel('3DM Threshold (mA)')
# g.axes.ravel()[1].set_ylabel('')
# plt.legend(loc='lower right')

#%% threshold coeff of variation intrafascicle and inter
sns.set(font_scale=1.75, style='whitegrid')
vardat = newdefdat.query('sim==3 and contact!="anodic"')
grouped = vardat.groupby(['sample', 'fiber_diam', 'sim', 'type', 'inner', 'deformation'])
analysis = grouped.agg({'threshold': [variation, np.mean]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)

plt.figure()
g = sns.catplot(
    data=analysis,
    y='threshold_variation',
    col='deformation',
    hue='fiber_diam',
    x='type',
    kind='box',
    palette='RdPu',
    facet_kws={'sharey': False},
)
for ax in g.axes.ravel():
    ax.set_yscale('log')
plt.suptitle('intrafascicle')
g.set_ylabels('Threshold Coeff. of Var.')
plt.subplots_adjust(top=0.85)
# plt.gca().set_aspect(20)
# plt.xlabel('')

# now do variance between fascicle mean thresholds
grouped = analysis.groupby(['sample', 'fiber_diam', 'type', 'deformation'])
analysis = grouped.agg({'threshold_mean': [variation]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)

plt.figure()
g = sns.catplot(
    data=analysis,
    y='threshold_mean_variation',
    col='deformation',
    hue='fiber_diam',
    x='type',
    kind='box',
    palette='RdPu',
    facet_kws={'sharey': False},
)
plt.suptitle('interfascicle')
g.set_ylabels('Threshold Coeff. of Var.')
plt.subplots_adjust(top=0.85)

# threshold variances for whole nerve
sns.set(font_scale=1.75, style='whitegrid')
vardat = newdefdat.query('sim==3 and contact!="anodic"')
grouped = vardat.groupby(['sample', 'fiber_diam', 'sim', 'type', 'deformation'])
analysis = grouped.agg({'threshold': [np.var, np.mean, variation]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.reset_index(inplace=True)
analysis.dropna(inplace=True)


plt.figure()
g = sns.catplot(
    data=analysis,
    y='threshold_variation',
    col='deformation',
    hue='fiber_diam',
    x='type',
    kind='box',
    palette='RdPu',
    facet_kws={'sharey': False},
)
plt.suptitle('intra-sample')
g.set_ylabels('Threshold Coeff. of Var.')
plt.subplots_adjust(top=0.85)

#%% peri-thk with 2DEM3DM
for stringdat in ['Undeformed', '2D-3D', '3D-3D']:
    subdat = newdefdat.query(f"deformation=='{stringdat}'")
    sns.set(font_scale=2)
    sns.set_style('whitegrid')
    nsimdata = subdat.query('fiber_diam in [3,13]')
    actlocdata = nsimdata.query('type=="3DM"')
    actlocdata['type'] = '3DM*'
    actlocdata['peri_thk'] = actlocdata['peri_thk_act_site']
    nsimdata = pd.concat([nsimdata, actlocdata])
    corrs = (
        nsimdata.groupby(['sample', 'fiber_diam', 'nerve_label', 'type'])['threshold', 'peri_thk'].corr().iloc[0::2, -1]
    )
    corrs = corrs.reset_index().rename(columns={'threshold': 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    means = corrs.groupby(['type', 'fiber_diam']).agg(np.mean)
    # plot one one line out to max of 3DM
    g = sns.lmplot(
        data=nsimdata.rename(columns={"nerve_label": "Sample"}),
        y='threshold',
        x='peri_thk',
        hue='Sample',
        palette=defpal,
        col='type',
        scatter_kws={'linewidths': 1, 'edgecolor': 'k', 's': 50},
        row='fiber_diam',
        facet_kws={'margin_titles': True, 'sharey': 'row'},
    )
    g.axes.flat[3].set_xlabel('Perineurium Thickness (μm)')
    g.axes.flat[4].set_xlabel('Perineurium Thickness (μm)')
    g.axes.flat[5].set_xlabel('Perineurium Thickness (μm)')
    g.axes.flat[0].set_ylabel('Threshold (mA)')
    g.axes.flat[3].set_ylabel('Threshold (mA)')

    for ax, r, letter in zip(
        g.axes.ravel(order='F'), means['peri_thk'], ['A', '', 'B', '', 'C', '']
    ):  # TODO: these are not ordering correctly
        nl = "\n"
        ax.text(0, 0.9, fr'{letter}{nl}$\bar{{r}}$={r:.2f}', transform=ax.transAxes)
        ax.set_ylim([0, None])
    # plt.xlabel('3DM Threshold (mA)')
    g.set_titles(col_template='{col_name}', row_template='Fiber Diameter = {row_name}')
    # g.fig.set_size_inches(10, 6)
    plt.suptitle("Deformation: " + stringdat, y=1.02)
#%%
defdr = newdefdat.copy()
defdr['percent_activated'] = 0
defdr = defdr.rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
for i, row in defdr.iterrows():
    thisdat = defdr.query(
        'modeltype == @row.modeltype and samplenum == @row.samplenum and fiber_diam == @row.fiber_diam and sim == @row.sim and deformation == @row.deformation'
    )
    # percent is number of thresholds less than or equal to this threshold divided by total number of thresholds
    defdr.loc[i, 'percent_activated'] = len(thisdat.query('threshold <= @row.threshold')) / len(thisdat)
defdr.sort_values('modeltype', inplace=True)
#%% activation order
for stringdat in ['Undeformed', '2D-3D', '3D-3D']:
    newdefdr = defdr.copy().query(f"deformation=='{stringdat}'")
    sns.set(font_scale=2)
    sns.set_style('whitegrid')
    newdefdr['percent_activated2d'] = np.nan
    # go through every row and for each fiber find the 2D activation percent
    for i, row in newdefdr.iterrows():
        # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
        thisdat = newdefdr.query(
            'modeltype == "2DEM" and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim and master_fiber_index == @row.master_fiber_index and contact == "cathodic"'
        )
        assert len(thisdat) == 1
        val = thisdat.percent_activated.values[0]
        assert not val is np.nan
        newdefdr.loc[i, 'percent_activated2d'] = val
    sns.set(font_scale=2.25, style='whitegrid')
    plt.figure()
    g = sns.catplot(
        kind='swarm',
        row='fiber_diam',
        data=newdefdr.query("fiber_diam in [3,13] and contact != 'anodic'"),
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
        s=25,
    )
    plt.subplots_adjust(top=0.87)
    plt.suptitle(stringdat, x=0.37)
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
        sm, ax=g.axes.ravel().tolist(), aspect=10, shrink=0.8, label='Proportion 2DEM fibers activated'
    ).ax.yaxis.set_ticks_position('left')

    # g._legend.set_title('')
#%% dose-response deformed
sns.set(font_scale=2)
sns.set_style('whitegrid')
# could shade in min and max to show response range?

plt.figure()
g = sns.relplot(
    kind='line',
    row='fiber_diam',
    data=defdr.query("fiber_diam in [3,13] and contact != 'anodic'"),
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
#%%
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
g = sns.relplot(data=allpemean, kind='line', col='deformation', style='level', y='pe', x='fiber_diam', markers=True)
g.legend.set_title('')
g.set_xlabels('Fiber Diameter (μm)')
g.set_titles(col_template='Deformation: {col_name}')
g.set_ylabels("Percent Difference")
plt.xticks([3, 5, 7, 9, 11, 13])
plt.ylim(ylimpe)

#%% organization compare type
for stringdat in ['Undeformed', '3D-3D']:
    # remove all samples which dont end with a 2
    sns.reset_orig()
    sns.set(font_scale=1.75, style='whitegrid')
    mpl.rcParams['figure.dpi'] = 400
    threshdat = newdefdat[~newdefdat['sample'].astype(str).str.endswith('0')]
    # subtract 1 from sample if type is 3DM
    threshdat['sample'] = threshdat.apply(lambda x: x['sample'] - 1 if x.type == "3DM" else x['sample'], axis=1)
    data = threshdat.query(f'nsim in [0,5] and sim == 3 and deformation=="{stringdat}"')
    # apply minmax normalization within sample and nsim
    g = sns.FacetGrid(
        data.rename(columns={"threshold": "Threshold (mA)", "nerve_label": "Sample"}),
        col="Sample",
        row='fiber_diam',
        sharey=False,
        margin_titles=True,
    )
    g.map_dataframe(
        sns.lineplot,
        x='type',
        y='Threshold (mA)',
        units='master_fiber_index',
        estimator=None,
        color='k',
        linewidth=1,
        alpha=0.25,
    )
    g.map_dataframe(sns.stripplot, jitter=False, linewidth=1, x='type', y='Threshold (mA)', palette='binary')

    plt.subplots_adjust(top=0.9)
    g.fig.set_size_inches(10, 6)
    g.set_titles(col_template='{col_name}', row_template='Fiber Diameter: {row_name}')
    g.set_ylabels('Threshold (mA)')
    plt.suptitle(f'thresholds compared between 2DEM (cathodic) and 3DM for {stringdat}')
#%% organization compare nsim
# remove all samples which dont end with a 2
for stringdat in ['Undeformed', '3D-3D']:
    sns.reset_orig()
    sns.set(font_scale=1.75, style='whitegrid')
    mpl.rcParams['figure.dpi'] = 400
    threshdat = newdefdat[~newdefdat['sample'].astype(str).str.endswith('0')]
    data = threshdat.query(f'nsim in [0,5] and sim == 3 and deformation=="{stringdat}"')
    # apply minmax normalization within sample and nsim
    data['Threshold (mA)'] = data.groupby(['sample', 'nsim', 'type'])['threshold'].transform(
        lambda x: (x - x.min()) / (x.max() - x.min())
    )
    # set nsim as categorical
    data['Fiber Diameter'] = data['nsim'].astype('category')
    # nsim 0 is three and nsim 5 is thirteen
    data['Fiber Diameter'].cat.rename_categories({0: '3', 5: '13'}, inplace=True)
    g = sns.FacetGrid(data, col="nerve_label", row='type', sharey=False, margin_titles=True)
    g.map_dataframe(
        sns.lineplot,
        x='Fiber Diameter',
        y='Threshold (mA)',
        units='master_fiber_index',
        estimator=None,
        color='k',
        alpha=0.25,
    )
    g.map_dataframe(sns.stripplot, linewidth=1, x='Fiber Diameter', y='Threshold (mA)', palette='binary', jitter=False)
    # plt.subplots_adjust(top=0.9)
    # plt.suptitle('min-max normalized thresholds compared between 3 um and 13 um thresholds (cath 2DEMEM)')
    g.set_titles(col_template='{col_name}', row_template='{row_name}')
    plt.savefig('matchsim.png', dpi=400)
    g.set_ylabels('Threshold (mA)')
    plt.gcf().set_size_inches([10, 6])
    plt.suptitle(f'thresholds compared between fiber diameters {stringdat}')
    plt.subplots_adjust(top=0.9)
#%% organizatino compare nsim
for stringdat in ['Undeformed', '3D-3D']:
    sns.set(font_scale=1.75)
    sns.set_style('whitegrid')
    threshdat = newdefdat[~newdefdat['sample'].astype(str).str.endswith('0')]

    scores = []
    for nerve in pd.unique(threshdat['nerve_label']):
        for n in [3, 13]:
            shortdat = threshdat.query(f'nerve_label=="{nerve}" and fiber_diam=={n} and deformation=="{stringdat}"')
            data2d = shortdat.query('type=="2DEM"').sort_values('threshold').master_fiber_index
            data3d = shortdat.query('type=="3DM"').sort_values('threshold').master_fiber_index
            rc = compute_reorder_cost(list(data2d), list(data3d))
            scores.append({'sample': nerve, 'fiber_diam': n, 'score2d3d': rc})
    plt.figure()
    scoredat = pd.DataFrame(scores)
    scoredat['fiber_diam'] = pd.Categorical(scoredat['fiber_diam'].astype(int), categories=[3, 13], ordered=True)
    ax = sns.boxplot(data=scoredat, x='score2d3d', y='fiber_diam', boxprops={'facecolor': 'white'})
    ax.set_ylabel('Fiber Diameter (μm)')
    plt.xlabel('RC')
    plt.gcf().set_size_inches([3, 5])
    plt.title(f'{stringdat}')
    plt.xlim([0, 0.5])
#%% organizatino compare type
for stringdat in ['Undeformed', '3D-3D']:
    sns.set(font_scale=1.75)
    sns.set_style('whitegrid')
    threshdat = newdefdat[~newdefdat['sample'].astype(str).str.endswith('0')]
    scores = []
    for nerve in pd.unique(threshdat['nerve_label']):
        for model in ['2DEM', '3DM']:
            shortdat = threshdat.query(f'nerve_label=="{nerve}" and type=="{model}" and deformation=="{stringdat}"')
            datasmol = shortdat.query('nsim==0').sort_values('threshold').master_fiber_index
            databeeg = shortdat.query('nsim==5').sort_values('threshold').master_fiber_index
            rc = compute_reorder_cost(list(datasmol), list(databeeg))
            scores.append({'sample': nerve, 'type': model, 'scoresmolbeeg': rc})
    plt.figure()
    scoredat = pd.DataFrame(scores)
    ax = sns.boxplot(data=scoredat, x='scoresmolbeeg', y='type', boxprops={'facecolor': 'white'})
    plt.xlabel('RC')
    plt.ylabel('')
    plt.gcf().set_size_inches([3, 5])
    plt.title(f'{stringdat}')
    plt.xlim([0, 0.5])
#%% tortuosity lineplot
for stringdat in ['Undeformed', '3D-3D']:
    sns.set(font_scale=2)
    sns.set_style('whitegrid')
    nsimdata = newdefdat.query(f'fiber_diam in [3,13] and type=="3DM" and deformation=="{stringdat}"')
    corrs = (
        nsimdata.groupby(['sample', 'fiber_diam', 'nerve_label', 'type'])['threshold', 'tortuosity']
        .corr()
        .iloc[0::2, -1]
    )
    corrs = corrs.reset_index().rename(columns={'threshold': 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
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
        ax.text(0.5, 0.9, f'\nmean r={r:.2f}', transform=ax.transAxes)
    g.axes.ravel()[0].set_ylabel('Threshold (mA)')
    # plt.xlabel('3DM Threshold (mA)')
    g.set_titles(col_template='Fiber Diameter = {col_name} μm')
    plt.suptitle(stringdat)
    plt.subplots_adjust(top=0.8)
    plt.show()

#%%begin 3d vs 3D code
#%% threshold violinplot
# generate boxplot of 2DEM and 3DM thresholds
sns.set(font_scale=2)
sns.set_style('whitegrid')
g = sns.catplot(
    kind='violin',
    col='nerve_label',
    row='fiber_diam',
    data=defdefcomp.query('fiber_diam in [3,13]'),
    x='deformed',
    y='threshold',
    sharey=False,
    palette=pal2d3d,
    margin_titles=True,
    split=False,
)
g.set_titles(col_template="{col_name}", row_template="D: {row_name} μm")
g.set_ylabels("Threshold (mA)")
plt.show()
#%% sample thresholds with unity line
defdefmatched = datamatch(
    defdefcomp.query('deformed == True'), defdefcomp.query('deformed == False'), 'threshold'
).drop(columns='deformed')
defdefmatched.rename(columns={'threshold3d': 'threshold', 'threshold': 'thresholddef'}, inplace=True)

sns.set(font_scale=1.75)
sns.set_style('whitegrid')
# usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
usedata = defdefmatched
nsimdata = usedata.query('fiber_diam in [3,13]')
g = sns.relplot(
    data=nsimdata,
    kind='scatter',
    col='fiber_diam',
    x='threshold',
    y='thresholddef',
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
    r, p = pearsonr(rdata.threshold, rdata.thresholddef)
    # perc = sum(rdata.threshold > rdata.thresholddef) / len(rdata.threshold)
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} μm: {r**2:.2f}')
    ax.set_title(f'{diam} μm')
    ax.plot([0, rdata.threshold.max()], [0, rdata.threshold.max()], '--k', linewidth=2, label='1:1 line')
    ax.set_xlabel('Undeformed Threshold (mA)')
    ax.set_yticks(ax.get_xticks())
    ax.set_aspect('equal', 'box')
    ax.set_xlim([0, None])
    ax.set_ylim([0, ax.get_xlim()[1]])
g.axes.ravel()[0].set_ylabel('Deformed Threshold (mA)')
g.axes.ravel()[1].set_ylabel('')
plt.legend(loc='lower right')
g.fig.set_size_inches([9.5, 5])
#%%
defdefmatched = datamatch(
    defdefcomp.query('deformed == True'), defdefcomp.query('deformed == False'), 'threshold'
).drop(columns='deformed')
defdefmatched.rename(columns={'threshold3d': 'threshold', 'threshold': 'thresholddef'}, inplace=True)

sns.set(font_scale=1.75)
sns.set_style('whitegrid')
# usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
usedata = defdefmatched
nsimdata = usedata.query('fiber_diam in [3,13]')
g = sns.lmplot(
    data=nsimdata,
    row='fiber_diam',
    x='threshold',
    y='thresholddef',
    col='nerve_label',
    # color='white',
    # s=20,
    palette=defpal,
    facet_kws={'sharex': False, 'sharey': False, "margin_titles": True},
    hue='nerve_label',
)
# g.map(sns.regplot, 'threshold', 'threshold3d', scatter=False, color='k', label='linear fit')
# plot one one line out to max of 3DM


#%% tortuosity lineplot
sns.set(font_scale=2)
sns.set_style('whitegrid')
nsimdata = defdefcomp.query(f'fiber_diam in [3,13]')
corrs = (
    nsimdata.groupby(['sample', 'fiber_diam', 'nerve_label', 'deformed'])['threshold', 'tortuosity']
    .corr()
    .iloc[0::2, -1]
)
corrs = corrs.reset_index().rename(columns={'threshold': 'correlation'})
corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
means = corrs.groupby(['deformed', 'fiber_diam']).agg(np.mean)
# plot one one line out to max of 3DM
g = sns.lmplot(
    data=nsimdata.rename(columns={"nerve_label": "Sample"}),
    y='threshold',
    x='tortuosity',
    hue='Sample',
    palette='colorblind',
    col='fiber_diam',
    row='deformed',
    scatter_kws={'linewidths': 1, 'edgecolor': 'k', 's': 50},
    facet_kws={'margin_titles': True, 'sharey': False},
)
# for ax, r in zip(g.axes.ravel(), means['tortuosity']):
# ax.text(0.5, 0.9, f'\nmean r={r:.2f}', transform=ax.transAxes)
g.set_ylabels('Threshold (mA)')
# plt.xlabel('3DM Threshold (mA)')
g.set_titles(col_template='Fiber Diameter = {col_name} μm')
plt.subplots_adjust(top=0.8)
plt.show()
#%% dose-response deformed
sns.set(font_scale=2)
sns.set_style('whitegrid')
# could shade in min and max to show response range?
defdefdr = defdr.query("modeltype=='3DM' and deformation!='3D-3D'")
defdefdr['deformed'] = defdefdr['deformation'] == "2D-3D"

plt.figure()
g = sns.relplot(
    kind='line',
    row='fiber_diam',
    data=defdefdr.query("fiber_diam in [3,13]"),
    y='percent_activated',
    x='threshold',
    hue='nerve_label',
    col="nerve_label",
    palette=defpal,
    estimator=None,
    linewidth=3,
    facet_kws={'sharex': False, 'margin_titles': True},
    **{'style': 'deformed'},
)
handles, labs = plt.gca().get_legend_handles_labels()
plt.legend(title='', handles=handles[1:], labels=labs[1:], bbox_to_anchor=[1, 1])

g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion fibers\nactivated (3 μm)')
g.axes[1][0].set_ylabel('Proportion fibers\nactivated (13 μm)')
g._legend.set_title('')
g.set_xlabels('Threshold (mA)')
g.set_titles(col_template="Deformation: {col_name}", row_template='')
g._legend.texts[0].set_text('')
g._legend.texts[5].set_text('')
#%%
newdefdef = newdefdat.query("type=='3DM' and deformation!='3D-3D'")
newdefdat['deformed'] = newdefdat['deformation'] == "2D-3D"

sns.set(font_scale=1.25)
sns.set_style('whitegrid')
plt.figure()
levels = {
    'onset': 10,
    'saturation': 90,
}
grouped = newdefdef.groupby(['sample', 'fiber_diam', 'sim', 'nerve_label', 'model', 'nsim', 'deformed'])
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
    id_vars=['sample', 'fiber_diam', 'sim', 'nerve_label', 'model', 'nsim', 'deformed'],
    value_vars=['onset', 'saturation'],
    var_name='level',
    value_name='threshold',
)

# set up facetgrid with nsim as row and level as columns
compiled_data.reset_index(inplace=True)
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
g.map_dataframe(sns.swarmplot, x='deformed', y='Threshold (mA)', color='black')
g.map_dataframe(
    sns.lineplot,
    x='deformed',
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
plt.suptitle('3D undeformed vs deformed')

compiled_data.dropna(inplace=True)

# calculate percent error for sample onset and saturation as well as population onset and saturation
pes = []
for level in ['onset', 'saturation']:
    for sam in compiled_data.nerve_label.unique():
        for fiber_diam in compiled_data.fiber_diam.unique():
            a = compiled_data.query(
                f"nerve_label == '{sam}' and level == '{level}' and deformed == False and fiber_diam == {fiber_diam}"
            )['threshold'].values
            b = compiled_data.query(
                f"nerve_label == '{sam}' and level == '{level}' and deformed == True and fiber_diam == {fiber_diam}"
            )['threshold'].values
            assert len(a) == len(b) == 1
            pe_res = pe(a[0], b[0])
            pes.append({'level': level, 'nerve_label': sam, 'fiber_diam': fiber_diam, 'pe': pe_res})
pes = pd.DataFrame(pes)
peses.append(pes)
# plot percent error for each sample DR
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
plt.figure()
sns.lineplot(data=pes, style='level', y='pe', x='fiber_diam', markers=True, hue='nerve_label')
plt.legend(bbox_to_anchor=[1, 1])
plt.xlabel('Fiber Diameter (μm)')
plt.ylabel('Percent Error')
# plt.title(f'Percent Error for {stringdat} Onset and Saturation')
plt.ylim([0, None])
ylimall = plt.gca().get_ylim()
# now calculate percent error for population onset and saturation
pemean = []
for level in ['onset', 'saturation']:
    for fiber_diam in compiled_data.fiber_diam.unique():
        a = compiled_data.query(f"level == '{level}' and deformed == False and fiber_diam == {fiber_diam}")[
            'threshold'
        ].values
        b = compiled_data.query(f"level == '{level}' and deformed == True and fiber_diam == {fiber_diam}")[
            'threshold'
        ].values
        pe_res = pe(np.mean(a), np.mean(b))
        pemean.append({'level': level, 'fiber_diam': fiber_diam, 'pe': pe_res})

pemean = pd.DataFrame(pemean)
pemeans.append(pemean)

# plot percent error for population
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
plt.figure()
sns.lineplot(data=pemean, style='level', y='pe', x='fiber_diam', markers=True)
plt.legend(title='')
plt.xlabel('Fiber Diameter (μm)')
plt.ylabel('Percent Error')
plt.title('Percent Diff 3D undeformed vs deformed')
plt.ylim(ylimall)

#%%
defdefdr = defdr.query("modeltype=='3DM' and deformation!='3D-3D'")
defdefdr['deformed'] = defdefdr['deformation'] == "2D-3D"
# mean across inners
defdefdr = defdefdr.groupby(['nerve_label', 'fiber_diam', 'sim', 'deformed', 'inner']).mean().reset_index()

defdefdr['percent_activatedUD'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for i, row in defdefdr.iterrows():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = defdefdr.query(
        'deformed == False and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim and inner == @row.inner'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert not val is np.nan
    defdefdr.loc[i, 'percent_activatedUD'] = val
#%%
sns.set(font_scale=2, style='whitegrid')
plt.figure()
g = sns.catplot(
    kind='swarm',
    row='fiber_diam',
    data=defdefdr.query("fiber_diam in [3,13]"),
    y='percent_activated',
    x='deformed',
    units='nerve_label',
    col='nerve_label',
    palette='plasma',
    hue='percent_activatedUD',
    # hue='inner',
    estimator=None,
    linewidth=0,
    facet_kws={'margin_titles': True},
    s=100,
)

g.set_titles(row_template='', col_template='{col_name}')
g.axes[0][0].set_xlabel('')
g.axes[0][0].set_ylabel('Proportion of fibers activated\nFiber diameter: 3 μm')
g.axes[1][0].set_ylabel('Proportion of fibers activated\nFiber diameter: 13 μm')
g.legend.set_title("Percent Activated Undeformed")
#%% organization compare type with recruitment order
# remove all samples which dont end with a 2
sns.reset_orig()
sns.set(font_scale=1.25, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
# threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]
# # subtract 1 from sample if type is 3DM
# threshdat['sample'] = threshdat.apply(lambda x: x['sample'] - 1 if x.type == "3DM" else x['sample'], axis=1)
data = defdefdr.query("fiber_diam in [3,13]")

# apply minmax normalization within sample and nsim
g = sns.FacetGrid(
    data,
    col="nerve_label",
    row='fiber_diam',
    sharey=False,
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x='deformed',
    y='percent_activated',
    units='master_fiber_index',
    estimator=None,
    color='k',
    linewidth=1,
    alpha=0.25,
)
g.map_dataframe(
    sns.stripplot,
    jitter=False,
    linewidth=1,
    x='deformed',
    y='percent_activated',
    palette='plasma',
    hue='percent_activatedUD',
)

plt.subplots_adjust(top=0.9)
g.fig.set_size_inches(10, 6)
g.set_titles(col_template='{col_name}', row_template='')
# g.set_xlabels('')
plt.legend(title="Percent Activated Undeformed", bbox_to_anchor=[1, 1])
g.axes.ravel()[0].set_ylabel('Percent\nFiber Diam: 3μm')
g.axes.ravel()[6].set_ylabel('Percent\nFiber Diam: 13μm')

# plt.suptitle('thresholds compared between 2DEM (cathodic) and 3DM')
g = sns.FacetGrid(
    data,
    col="nerve_label",
    row='fiber_diam',
    sharey=False,
    margin_titles=True,
)
g.map_dataframe(
    sns.lineplot,
    x='deformed',
    y='threshold',
    units='master_fiber_index',
    estimator=None,
    color='k',
    linewidth=1,
    alpha=0.25,
)
g.map_dataframe(
    sns.stripplot, jitter=False, linewidth=1, x='deformed', y='threshold', palette='plasma', hue='percent_activatedUD'
)

plt.subplots_adjust(top=0.9)
g.fig.set_size_inches(10, 6)
g.set_titles(col_template='{col_name}', row_template='')
# g.set_xlabels('')
plt.legend(title="Percent Activated 2DEM", bbox_to_anchor=[1, 1])
g.axes.ravel()[0].set_ylabel('Threshold\nFiber Diam: 3μm')
g.axes.ravel()[6].set_ylabel('Threshold\nFiber Diam: 13μm')
# plt.suptitle('thresholds compared between 2DEM (cathodic) and 3DM')
sns.set(font_scale=1.75)
sns.set_style('whitegrid')

#%%
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

#%% TODO with deformed
# generate boxplot of 2DEM and 3DM thresholds
# sns.set(font_scale=1.75)
sns.set_style('whitegrid')
g = sns.catplot(
    kind='violin',
    col='fiber_diam',
    row='nerve_label',
    data=threshload,
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
# g.fig.set_size_inches(5, 5, forward=True)
plt.subplots_adjust(top=0.9, wspace=0.4)
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# plt.suptitle('2DEM vs 3DM Thresholds for all samples')
g.axes.ravel()[0].set_title('3 μm')
g.axes.ravel()[1].set_title('13 μm')
plt.show()
# calculate mean, std, and median for each sample
outdata = matched.groupby(['nerve_label', 'fiber_diam']).agg(
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
plt.title('Mean Threshold Percent Difference')
plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0))
plt.show()
sns.lineplot(data=outdata, x='fiber_diam', y='std_diff', hue='nerve_label', palette='colorblind', marker='o')
plt.ylabel('Percent Difference')
plt.xlabel('Fiber Diam')
plt.title('Standard Deviation Percent Difference')
plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0))
plt.show()
sns.lineplot(data=outdata, x='fiber_diam', y='median_diff', hue='nerve_label', palette='colorblind', marker='o')
plt.ylabel('Percent Difference')
plt.xlabel('Fiber Diam')
plt.title('Median Threshold Percent Difference')
plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0))
plt.show()
