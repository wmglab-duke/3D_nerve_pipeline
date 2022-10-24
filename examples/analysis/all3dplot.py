# -*- coding: utf-8 -*-
import json
import os

import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import pearsonr, variation

os.chdir('../../')
import sys

sys.path.pop(-2)
from src.core.plotter import datamatch, datamatchlist

mpl.rcParams['figure.dpi'] = 400


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


#%% Setup
sns.set_style('whitegrid')
# base data
threshload = pd.read_csv('thresh_unmatched_sim3.csv').query('sim==3')
threshload = addpwfd(threshload, '3')
threshload['fiber_diam'] = threshload['fiber_diam'].astype(int)
# 3d and 2d trhresholds as different col same row
matched = datamatch(threshload.query('type=="2D"'), threshload.query('type=="3D"'), 'threshold').drop(columns='type')
# add new EMsample column to matched dataframe, which is the nerve label plus the first letter of the contact type capitalized
matched['EMsample'] = matched['nerve_label'] + matched['contact'].str[0].str.upper()
# for matched boxplots
repeated = matched.melt(
    id_vars=[x for x in matched.columns if 'threshold' not in x],
    value_vars=['threshold3d', 'threshold'],
    var_name='type',
    value_name='threshold',
)
repeated['type'] = repeated['type'].replace({'threshold': '2D', 'threshold3d': '3D'})
repeated['type'] = pd.Categorical(repeated['type'], categories=['2D', '3D'], ordered=True)

multimatched = datamatchlist(
    threshload.query('type=="2D"'), threshload.query('type=="3D"'), ['threshold', 'activation_zpos']
).drop(columns='type')

sys.exit()

#%% all thresholds with unity line
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
sns.scatterplot(data=matched.sort_values('fiber_diam'), x='threshold', y='threshold3d', hue='fiber_diam', s=20)
# plot one one line out to max of 3d
plt.plot([0, matched.threshold.max()], [0, matched.threshold.max()], 'r', linewidth=2)
plt.ylabel('3D Threshold (mA, log scale)')
plt.xlabel('2D Threshold (mA, log scale)')
# set legend title
plt.legend(title='Fiber Diameter (μm)')
plt.xscale('log')
plt.yscale('log')
# make axes squre
plt.gca().set_aspect('equal', adjustable='box')
# calculate correlation between 2d and 3d and percentage 2d greater than 3d, add to title
r, p = pearsonr(matched.threshold, matched.threshold3d)
perc = sum(matched.threshold > matched.threshold3d) / len(matched.threshold)
plt.title(
    '2D vs 3D Thresholds for all samples and fiber diameters\n'
    f'Correlation between 2D and 3D thresholds: r={r:.2f}, p={p:.2f}\n'
    f'percentage of 3D threshold lower than 2D counterpart: {perc:.2f}'
)
sns.move_legend(plt.gca(), "upper left", bbox_to_anchor=(1, 1))
plt.show()
#%% sample thresholds with unity line
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
for diam in [3, 13]:
    nsimdata = matched.query(f'fiber_diam=={diam}')
    sns.scatterplot(data=nsimdata, x='threshold', y='threshold3d', hue='nerve_label', s=20, palette='colorblind')
    # plot one one line out to max of 3d
    plt.plot([0, nsimdata.threshold.max()], [0, nsimdata.threshold.max()], 'r', linewidth=2)
    plt.ylabel('3D Threshold (mA)')
    plt.xlabel('2D Threshold (mA)')
    # set legend title
    plt.legend(title='Sample')
    # plt.xscale('log')
    # plt.yscale('log')
    # # make axes squre
    plt.gca().set_aspect('equal', adjustable='box')
    sns.move_legend(plt.gca(), "upper left", bbox_to_anchor=(1, 1))
    r, p = pearsonr(nsimdata.threshold, nsimdata.threshold3d)
    perc = sum(nsimdata.threshold > nsimdata.threshold3d) / len(nsimdata.threshold)
    # plt.title(
    #     f'2D vs 3D Thresholds for all samples and {diam} micron fibers\n'
    #     f'Correlation between 2D and 3D thresholds: r={r:.2f}, p={p:.2f}\n'
    #     f'percentage of 3D threshold lower than 2D counterpart: {perc:.2f}'
    # )
    # add correlation to plot
    plt.text(1.08, 0.2, f'r={r:.2f}', transform=plt.gca().transAxes)
    plt.title(f'Fiber Diameter: {diam} μm')
    plt.show()
#%% sample thresholds with unity line
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
usedata = matched = datamatch(usedata.query('type=="2D"'), usedata.query('type=="3D"'), 'threshold').drop(
    columns='type'
)
nsimdata = usedata.query('fiber_diam in [3,13]')
g = sns.relplot(
    data=nsimdata,
    kind='scatter',
    row='fiber_diam',
    x='threshold',
    y='threshold3d',
    hue='nerve_label',
    s=20,
    palette='colorblind',
    facet_kws={'sharex': False, 'sharey': False},
)
# plot one one line out to max of 3d

# set legend title
# plt.legend(title='Sample')
# plt.xscale('log')
# plt.yscale('log')
# # make axes squre
# plt.gca().set_aspect('equal', adjustable='box')
# sns.move_legend(plt.gca(), "upper left", bbox_to_anchor=(1, 1))
for diam, pos, ax in zip([3, 13], (0.2, 0.8), g.axes.ravel()):
    rdata = nsimdata.query(f'fiber_diam=={diam}')
    r, p = pearsonr(rdata.threshold, rdata.threshold3d)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    # plt.title(
    #     f'2D vs 3D Thresholds for all samples and {diam} micron fibers\n'
    #     f'Correlation between 2D and 3D thresholds: r={r:.2f}, p={p:.2f}\n'
    #     f'percentage of 3D threshold lower than 2D counterpart: {perc:.2f}'
    # )
    # add correlation to plot
    ax.text(0, 1, f'r={r:.2f}', transform=ax.transAxes)
    ax.plot([0, rdata.threshold.max()], [0, rdata.threshold.max()], 'r', linewidth=2)
    ax.set_ylabel('3D Threshold (mA)')
    ax.set_xlabel('2D Threshold (mA)')
#%% threshold barplot
# generate boxplot of 2d and 3d thresholds
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
g = sns.catplot(
    kind='bar',
    col='fiber_diam',
    data=threshload.query("nsim in [0,5]"),
    x='type',
    y='threshold',
    sharey=False,
    errorbar='se',
    palette='colorblind',
)
g.fig.set_size_inches(4, 6, forward=True)
plt.subplots_adjust(top=0.9, wspace=0.4)
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
# plt.suptitle('2D vs 3D Thresholds for all samples')
g.axes.ravel()[0].set_title('3 μm')
g.axes.ravel()[1].set_title('13 μm')
plt.show()
#%% threshold barplot
# generate boxplot of 2d and 3d thresholds
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
g = sns.barplot(
    data=threshload.query("nsim in [0]"),
    x='type',
    y='threshold',
    errorbar='se',
    palette='binary',
)
plt.gca().set_aspect(2)
plt.subplots_adjust(top=0.9, wspace=0.4)
plt.ylabel('Threshold (mA)')
plt.show()

#%% organization compare nsim
# remove all samples which dont end with a 2
sns.reset_orig()
mpl.rcParams['figure.dpi'] = 400
threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]
data = threshdat.query('nsim in [0,5] and sim == 3')
# apply minmax normalization within sample and nsim
data['threshold'] = data.groupby(['sample', 'nsim', 'type'])['threshold'].transform(
    lambda x: (x - x.min()) / (x.max() - x.min())
)
# set nsim as categorical
data['nsim'] = data['nsim'].astype('category')
# nsim 0 is three and nsim 5 is thirteen
data['nsim'].cat.rename_categories({0: '3', 5: '13'}, inplace=True)
g = sns.FacetGrid(data, col="nerve_label", row='type', sharey=False, margin_titles=True)
g.map_dataframe(sns.boxplot, x='nsim', y='threshold', palette='colorblind', boxprops={'facecolor': "none"})
g.map_dataframe(sns.stripplot, linewidth=1, x='nsim', y='threshold', palette='colorblind', jitter=False)
g.map_dataframe(sns.lineplot, x='nsim', y='threshold', units='master_fiber_index', estimator=None, color='k')
plt.subplots_adjust(top=0.9)
plt.suptitle('min-max normalized thresholds compared between 3 um and 13 um thresholds (cath 2DEM)')
plt.savefig('matchsim.png', dpi=400)
#%% organization compare type
# remove all samples which dont end with a 2
sns.reset_orig()
sns.set(font_scale=1.25, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]
# subtract 1 from sample if type is 3D
threshdat['sample'] = threshdat.apply(lambda x: x['sample'] - 1 if x.type == "3D" else x['sample'], axis=1)
data = threshdat.query('nsim in [0,5] and sim == 3')
# apply minmax normalization within sample and nsim
g = sns.FacetGrid(
    data.rename(columns={"threshold": "Threshold (mA)", "nerve_label": "Sample", "fiber_diam": "Fiber Diameter (μm)"}),
    col="Sample",
    row='Fiber Diameter (μm)',
    sharey=False,
    margin_titles=True,
)
g.map_dataframe(sns.boxplot, x='type', y='Threshold (mA)', palette='colorblind', boxprops={'facecolor': "none"})
g.map_dataframe(sns.stripplot, jitter=False, linewidth=1, x='type', y='Threshold (mA)', palette='colorblind')
g.map_dataframe(
    sns.lineplot,
    x='type',
    y='Threshold (mA)',
    units='master_fiber_index',
    estimator=None,
    color='k',
    alpha=0.25,
    linewidth=1,
)
plt.subplots_adjust(top=0.9)
# g.fig.set_size_inches(12, 8)
# plt.suptitle('thresholds compared between 2D (cathodic) and 3D')
plt.savefig('matchsim.png', dpi=400)
#%% All thresholds compare
g = sns.catplot(
    data=repeated.query('nsim in [0,5]'),
    x='EMsample',
    y='threshold',
    hue='type',
    kind='box',
    palette='colorblind',
    col='fiber_diam',
    sharey=False,
)
g.set(ylim=(0, None))
plt.subplots_adjust(top=0.9)
plt.suptitle('2D vs 3D Thresholds for all samples')
#%% coeff of variation
grouped = repeated.groupby(['sample', 'fiber_diam', 'sim', 'type'])
# calculate coefficient of variation for thresholds on grouped data
analysis = grouped.agg({'threshold': [variation]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis = analysis.reset_index()
sns.catplot(
    data=analysis.rename(columns={"fiber_diam": "Fiber Diameter(μm)"}),
    x='type',
    y='threshold_variation',
    hue='Fiber Diameter(μm)',
    kind='box',
    palette='RdPu',
)
plt.ylabel('Threshold Coefficient of Variation')
# plt.title('Coefficient of Variation of Thresholds for 2D and 3D')
#%% threshold variances intrafascicle and inter
vardat = repeated.query('nsim in [0,5] and sim==3')
grouped = vardat.groupby(['EMsample', 'fiber_diam', 'sim', 'type', 'inner'])
analysis = grouped.agg({'threshold': [np.var, np.mean]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.reset_index(inplace=True)
g = sns.catplot(
    errorbar='se',
    data=analysis,
    sharey=False,
    col='fiber_diam',
    x='EMsample',
    y='threshold_var',
    hue='type',
    kind='bar',
    palette='colorblind',
)
# plt.suptitle('Variance values of thresholds within fascicles')
g.axes.ravel()[0].set_ylabel('Variance (mA^2)')
g = sns.catplot(
    errorbar='se',
    data=analysis,
    sharey=False,
    x='type',
    col='fiber_diam',
    y='threshold_var',
    kind='bar',
    palette='colorblind',
    aspect=0.5,
)
g.axes.ravel()[0].set_ylabel('Variance (mA^2)')
plt.subplots_adjust(top=0.9)
for ax, diam in zip(g.axes.ravel(), [3, 13]):
    ax.set_xlabel('')
    ax.set_title(f'{diam} μm')
# plt.suptitle('Variance values of thresholds within fascicles')
plt.figure()
g = sns.barplot(
    errorbar='se',
    data=analysis.query("fiber_diam==3"),
    y='threshold_var',
    x='type',
    palette='binary',
)
plt.ylabel('Variance (mA^2)')
plt.gca().set_aspect(100)
# now do variance between fascicle mean thresholds
grouped = analysis.groupby(['EMsample', 'fiber_diam', 'type'])
analysis = grouped.agg({'threshold_mean': [np.var]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.reset_index(inplace=True)
g = sns.catplot(
    data=analysis,
    sharey=False,
    col='fiber_diam',
    x='EMsample',
    y='threshold_mean_var',
    hue='type',
    kind='bar',
    palette='colorblind',
)
g.axes.ravel()[0].set_ylabel('Variance (mA^2)')
# plt.suptitle('Variance values of fascicle thresholds')
g = sns.catplot(
    data=analysis,
    sharey=False,
    x='type',
    col='fiber_diam',
    y='threshold_mean_var',
    kind='bar',
    palette='colorblind',
    aspect=0.5,
    errorbar=('se'),
)
for ax, diam in zip(g.axes.ravel(), [3, 13]):
    ax.set_xlabel('')
    ax.set_title(f'{diam} μm')
g.axes.ravel()[0].set_ylabel('Variance (mA^2)')
plt.subplots_adjust(top=0.9)
# plt.suptitle('Variance values of fascicle thresholds')
plt.figure()
g = sns.barplot(
    errorbar='se',
    data=analysis.query("fiber_diam==3"),
    y='threshold_mean_var',
    x='type',
    palette='binary',
)
plt.ylabel('Variance (mA^2)')
plt.gca().set_aspect(10)
#%% Dose-response info
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
plt.figure()
levels = {
    'onset': 0.01,
    'half': 0.5,
    'saturation': 1,
}
grouped = threshload.query('sim==3 and nsim in [0,5]').groupby(
    ['sample', 'fiber_diam', 'sim', 'type', 'nerve_label', 'model', 'nsim']
)
analysis = grouped.agg({'threshold': [np.amin, np.median, np.amax]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.rename(
    columns={'threshold_amin': 'onset', 'threshold_median': 'half', 'threshold_amax': 'saturation'}, inplace=True
)
analysis = analysis.reset_index()
# combine onset, saturation, and half into one column with identifier
compiled_data = analysis.melt(
    id_vars=['sample', 'fiber_diam', 'sim', 'type', 'nerve_label', 'model', 'nsim'],
    value_vars=['onset', 'half', 'saturation'],
    var_name='level',
    value_name='threshold',
)
# subtract 1 from index if type is 3D
for fiber_diam in pd.unique(compiled_data.fiber_diam):
    ax = sns.boxplot(
        data=compiled_data.query(f'fiber_diam=={int(fiber_diam)}'),
        x='level',
        y='threshold',
        hue='type',
        palette='colorblind',
        color='w',
        boxprops={'facecolor': "none"},
    )
    sns.swarmplot(
        data=compiled_data.query(f'fiber_diam=={int(fiber_diam)}'),
        x='level',
        y='threshold',
        hue='type',
        palette='colorblind',
        dodge=True,
    )
    ax.set_ylim([0, None])
    ax.set_ylabel('Threshold (mA)')
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[2:], labels[2:])
    plt.title(f'fiber_diam={int(fiber_diam)} μm')
    plt.show()
# TODO get this to be paired boxplots again an dline up datapoints
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
    compiled_data.rename(columns={'fiber_diam': "Fiber Diameter (μm)", 'threshold': "Threshold (mA)"}),
    col="level",
    row="Fiber Diameter (μm)",
    sharey=False,
    margin_titles=True,
)
g.map_dataframe(sns.boxplot, x='type', y='Threshold (mA)', palette='colorblind', boxprops={'facecolor': "none"})
g.map_dataframe(sns.swarmplot, x='type', y='Threshold (mA)', palette='colorblind')
g.map_dataframe(sns.lineplot, x='type', y='Threshold (mA)', units='units', estimator=None, color='k')
plt.subplots_adjust(top=0.9)
# plt.suptitle('Dose-response curves changes between model types (2D is cathodic leading only)')
# get datamatch for compiled_data
plt.figure()
from src.core.plotter import datamatch_agg

errordr = datamatch_agg(compiled_data.query('type=="2D"'), compiled_data.query('type=="3D"'), 'threshold').drop(
    columns='type'
)
errordr['pe'] = errordr.apply(lambda row: pe(row['threshold3d'], row['threshold']), axis=1)
ax = sns.barplot(data=errordr, hue='fiber_diam', x='level', y='pe', errorbar='se')
plt.title('Percent error between 2D and 3D dose-response (per sample)')
plt.ylabel('Threshold Percent error')

plt.figure()
# now calculate mean population threshold for 2D and 3D
popdat = errordr.groupby(['fiber_diam', 'level'])['threshold', 'threshold3d'].mean()
popdat.reset_index(inplace=True)
popdat['pe'] = popdat.apply(lambda row: pe(row['threshold3d'], row['threshold']), axis=1)
popdat['level'] = pd.Categorical(popdat['level'], categories=['onset', 'half', 'saturation'], ordered=True)
plt.ylim(ax.get_ylim())
sns.barplot(data=popdat, hue='fiber_diam', x='level', y='pe', errorbar='se')
plt.title('Percent error between 2D and 3D dose-response (population mean)')
plt.ylabel('Threshold Percent error')

#%% all dose-response
# TODO add percent activated to initial processing
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
drdat = threshload.copy()
drdat['percent_activated'] = 0
drdat = drdat.rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
for i, row in drdat.iterrows():
    thisdat = drdat.query(
        'modeltype == @row.modeltype and samplenum == @row.samplenum and fiber_diam == @row.fiber_diam and sim == @row.sim'
    )
    # percent is number of thresholds less than or equal to this threshold divided by total number of thresholds
    drdat.loc[i, 'percent_activated'] = len(thisdat.query('threshold <= @row.threshold')) / len(thisdat)
    # could shade in min and max to show response range
plt.figure()
drdat.sort_values('modeltype', inplace=True)
g = sns.relplot(
    data=drdat.query("fiber_diam in [3,13]"),
    facet_kws=dict(sharey=False),
    kind='scatter',
    x='percent_activated',
    row='fiber_diam',
    y='threshold',
    hue='modeltype',
    palette='colorblind',
)
for ax, diam in zip(g.axes.ravel(), [3, 13]):
    ax.set_xlabel('Proportion of fibers activated')
    ax.set_title(f'{diam} μm')
    ax.set_ylim([0, None])
    ax.set_ylabel('Threshold (mA)')
g.fig.set_size_inches([6, 6])
plt.legend(title='')
plt.figure()
g = sns.scatterplot(
    data=drdat.query("fiber_diam in [3,13]"),
    x='percent_activated',
    y='threshold',
    hue='modeltype',
    style='fiber_diam',
    palette='colorblind',
)
plt.yscale('log')
sns.move_legend(plt.gca(), "upper left", bbox_to_anchor=(1, 1))
plt.ylabel('Threshold (mA)')
plt.xlabel('Proportion of fibers activated')
plt.figure()
g = sns.lineplot(
    data=drdat.query("fiber_diam in [3,13]"),
    x='percent_activated',
    y='threshold',
    hue='modeltype',
    style='fiber_diam',
    palette='colorblind',
)
plt.yscale('log')
# sns.move_legend(plt.gca(), "upper left", bbox_to_anchor=(1, 1))
plt.ylabel('Threshold (mA)')
plt.xlabel('Proportion of fibers activated')
# plt.suptitle(f'Dose-Response for all samples fiber diam={fiber_diam} μm')
plt.figure()
g = sns.lineplot(
    data=drdat.query("fiber_diam in [3] and contact != 'anodic'"),  # stopped here
    x='percent_activated',
    y='threshold',
    units='nerve_label',
    hue='modeltype',
    palette='colorblind',
    estimator=None,
)
# plt.yscale('log')
sns.move_legend(plt.gca(), "upper left", bbox_to_anchor=(1, 1))
plt.ylabel('Threshold (mA)')
plt.xlabel('Proportion of fibers activated')
#%% strength-duration
# seaborn facet scatterplot with x as pulse width, y as threshold, and column as diameter
threshsd = addpwfd(pd.read_csv('thresh_unmatched_sim7.csv'), '7')
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
g = sns.catplot(
    kind='violin',
    data=threshsd,
    col="fiber_diam",
    col_wrap=2,
    sharey=False,
    x='pulse_width',
    y='threshold',
    palette='colorblind',
    hue='type',
    dodge=True,
    linewidth=1,
)
plt.subplots_adjust(top=0.85)
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
plt.suptitle('Threshold vs. Pulse Width')
#%%z-score
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
threshdat = threshload.copy()
threshdat['z_score'] = 0
threshdat['z_score_error'] = 0
threshdat = threshdat.rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
for i, row in threshdat.iterrows():
    thisdat = threshdat.query(
        'modeltype == @row.modeltype and samplenum == @row.samplenum and fiber_diam == @row.fiber_diam and sim == @row.sim'
    )
    # percent is number of thresholds less than or equal to this threshold divided by total number of thresholds
    threshdat.loc[i, 'z_score'] = (row.threshold - np.mean(thisdat.threshold)) / np.std(thisdat.threshold)
for i, row in threshdat.iterrows():
    # subtract this row z score from the z score of fiber_diam 0
    threshdat.loc[i, 'z_score_error'] = float(
        row.z_score
        - np.mean(
            threshdat.query(
                'master_fiber_index == @row.master_fiber_index and modeltype == @row.modeltype and samplenum == @row.samplenum and sim == @row.sim'
            ).z_score
        )
    )
sns.catplot(data=threshdat, x='modeltype', y='z_score_error', hue='fiber_diam', kind='box', palette='colorblind')
plt.title('Z Score Residual for 2D and 3D')
plt.ylabel('Z Score Residual from mean')
#%% Plot histogram of activation z positions
sns.set(font_scale=1.5, style='whitegrid')
g = sns.displot(data=threshload, y='activation_zpos', row='fiber_diam', hue='type', kind='kde', palette='colorblind')
g.fig.set_size_inches([3, 10])
g = sns.displot(
    data=threshload.query("type=='3D' and nsim in [0,5]"),
    y='activation_zpos',
    col='fiber_diam',
    hue='nerve_label',
    facet_kws={'sharex': False},
    kind='kde',
    palette='colorblind',
)
g = sns.displot(
    data=threshload.query("type=='2D' and nsim in [0,5]"),
    y='activation_zpos',
    col='fiber_diam',
    hue='nerve_label',
    facet_kws={'sharex': False},
    kind='kde',
    palette='colorblind',
)
g = sns.displot(
    data=threshload.query("nsim in [0]").rename(columns={'nerve_label': 'Sample'}),
    y='activation_zpos',
    col='type',
    hue='Sample',
    facet_kws={'sharex': False},
    kind='hist',
    kde=True,
    palette='colorblind',
)
g.axes.ravel()[0].set_ylabel('Z-position of activation (μm)')
# plt.suptitle('Activation z-position for 3um fibers')

g = sns.displot(
    data=threshload.query("nsim in [5]").rename(columns={'nerve_label': 'Sample'}),
    y='activation_zpos',
    col='type',
    hue='type',
    facet_kws={'sharex': False},
    kind='hist',
    kde=True,
    stat='probability',
    palette='binary',
)
g.axes.ravel()[0].set_ylabel('Z-position of activation (μm)')
# plt.suptitle('Activation z-position for 3um fibers')
g = sns.displot(
    data=threshload.query("nsim in [5]").rename(columns={'nerve_label': 'Sample'}),
    y='activation_zpos',
    # col='type',
    hue='type',
    facet_kws={'sharex': False},
    kind='hist',
    kde=True,
    stat='probability',
    palette='colorblind',
)
plt.ylim([19000, 31000])
g.axes.ravel()[0].set_ylabel('Z-position of activation (μm)')
# plt.suptitle('Activation z-position for 3um fibers')
g = sns.displot(
    data=threshload.query("nsim in [0,5]").rename(columns={'nerve_label': 'Sample'}),
    y='activation_zpos',
    col='type',
    row='fiber_diam',
    hue='Sample',
    facet_kws={'sharex': False},
    kind='kde',
    palette='colorblind',
    common_norm=False,
)
ylim = g.axes[0][0].get_ylim()
zloadmono = pd.read_csv('thresh_unmatched_sim10.csv')
zloadmono = addpwfd(zloadmono, '10')
zloadmono['fiber_diam'] = zloadmono['fiber_diam'].astype(int)
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
)
for ax in g.axes.ravel():
    plt.ylim(ylim)

#%% Correlation2d3d
sns.reset_orig()
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
mpl.rcParams['figure.dpi'] = 400
usedata = matched.rename(columns={'threshold': 'threshold2d'})
for comparison in [
    ['threshold2d', 'threshold3d'],
    ['threshold2d', 'peri_thk'],
    ['threshold3d', 'peri_thk'],
]:
    corrs = usedata.groupby(['sample', 'fiber_diam', 'contact', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
    corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    corrs['contact'] = pd.Categorical(corrs['contact'], categories=['cathodic', 'anodic'], ordered=True)
    plt.figure()
    corrs.rename(inplace=True, columns={"nerve_label": "Sample"})
    sns.scatterplot(
        data=corrs, x='fiber_diam', y='correlation', hue='Sample', s=100, palette='colorblind', style='contact'
    )
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
    plt.gca().set_ylim([0, 1])
    plt.xlabel('Fiber Diameter (μm)')
    # plt.title(f'Correlation between {comparison[0]} and {comparison[1]}')

#%% 3D correlations
sns.reset_orig()
mpl.rcParams['figure.dpi'] = 400
usedata = threshload.rename(columns={'threshold': 'threshold3d'})
for comparison in [
    ['threshold3d', 'tortuosity'],
    ['threshold3d', 'peri_thk'],
    ['threshold3d', 'peri_thk_act_site'],
    ['threshold3d', 'smallest_thk_under_cuff'],
    ['threshold3d', 'cuff_tortuosity'],
]:
    corrs = (
        usedata.query('type=="3D"').groupby(['sample', 'fiber_diam', 'nerve_label'])[comparison].corr().iloc[0::2, -1]
    )
    corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
    corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
    plt.figure()
    sns.scatterplot(data=corrs, x='fiber_diam', y='correlation', hue='nerve_label', s=100, palette='colorblind')
    ax = sns.lineplot(
        data=corrs, x='fiber_diam', y='correlation', hue='nerve_label', legend=False, palette='colorblind'
    )
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.gca().set_ylim([-1, 1])
    plt.title(f'Correlation between {comparison[0]} and {comparison[1]}')
    #%% 3D correlations monopolar
    sns.reset_orig()
    mpl.rcParams['figure.dpi'] = 400
    usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10').rename(columns={'threshold': 'threshold3d'})
    for comparison in [
        ['threshold3d', 'tortuosity'],
        ['threshold3d', 'peri_thk'],
        ['threshold3d', 'peri_thk_act_site'],
        # ['threshold3d', 'smallest_thk_under_cuff'],
        ['threshold3d', 'cuff_tortuosity'],
    ]:
        corrs = (
            usedata.query('type=="3D"')
            .groupby(['sample', 'fiber_diam', 'nerve_label'])[comparison]
            .corr()
            .iloc[0::2, -1]
        )
        corrs = corrs.reset_index().rename(columns={comparison[1]: 'correlation'})
        corrs['fiber_diam'] = pd.Categorical(corrs['fiber_diam'].astype(int), ordered=True)
        plt.figure()
        sns.scatterplot(data=corrs, x='fiber_diam', y='correlation', hue='nerve_label', s=100, palette='colorblind')
        ax = sns.lineplot(
            data=corrs, x='fiber_diam', y='correlation', hue='nerve_label', legend=False, palette='colorblind'
        )
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        plt.gca().set_ylim([-1, 1])
        plt.title(f'Correlation between {comparison[0]} and {comparison[1]} (MM)')
    #%%correlations monopolar
    sns.set(font_scale=1.25)
    sns.set_style('whitegrid')
    usedata = addpwfd(pd.read_csv('thresh_unmatched_sim10.csv'), '10')
    usedata = datamatch(usedata.query('type=="2D"'), usedata.query('type=="3D"'), 'threshold').drop(columns='type')
    # add new EMsample column to matched dataframe, which is the nerve label plus the first letter of the contact type capitalized
    usedata['EMsample'] = usedata['nerve_label'] + usedata['contact'].str[0].str.upper()
    usedata = usedata.rename(columns={'threshold': 'threshold2d'})
    for comparison in [
        ['threshold2d', 'threshold3d'],
        ['threshold2d', 'peri_thk'],
        ['threshold3d', 'peri_thk'],
    ]:  # TODO: make this into a function
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
        plt.gca().set_ylim([0, 1])
        # plt.title(f'Correlation between {comparison[0]} and {comparison[1]} (MM)')
#%%
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
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
plt.figure()
sns.histplot(data=threshload.query('type=="3D"'), x='tortuosity')
plt.figure()
g = sns.lmplot(
    data=threshload.query('type=="3D" and nsim in [0,5]').rename(columns={'nerve_label': 'Sample'}),
    facet_kws={'sharex': False, 'sharey': False},
    col='fiber_diam',
    hue='Sample',
    y='threshold',
    x='tortuosity',
    palette='colorblind',
)
g.axes.ravel()[0].set_ylabel('Threshold (mA)')
g.axes.ravel()[0].set_title('3 μm')
g.axes.ravel()[1].set_title('13 μm')
plt.figure()
g = sns.lmplot(
    data=threshload.query('type=="3D"and nsim==0').rename(columns={'nerve_label': 'Sample'}),
    facet_kws={'sharex': False, 'sharey': False},
    # col='fiber_diam',
    hue='Sample',
    y='threshold',
    x='tortuosity',
    palette='colorblind',
)
plt.ylabel('Threshold (mA)')
#%% Percent Error
sns.reset_orig()
mpl.rcParams['figure.dpi'] = 400
# sns.set(font_scale=1.25)
# apply pe to all rows of dataframe matched, with threshold3d as the correct value and threshold as the estimated value
multimatched['pe'] = multimatched.apply(lambda row: pe(row['threshold3d'], row['threshold']), axis=1)
# calculate difference between activation_zpos and activation_zpos3d
multimatched['zdiff'] = multimatched['activation_zpos'] - multimatched['activation_zpos3d']
multimatched['zdiff_abs'] = multimatched['zdiff'].abs()
sns.barplot(data=multimatched, x='fiber_diam', y='pe', errorbar='se')
plt.title('Threshold Percent Error by fiber diameter')
plt.ylabel('Percent Error')
# barplot of percent error by nerve label and fiber diameter
plt.figure()
sns.barplot(data=multimatched, x='nerve_label', y='pe', errorbar='se')
plt.title('Threshold Percent Error by sample')
plt.ylabel('Percent Error')
plt.figure()
sns.barplot(data=multimatched, x='nerve_label', y='pe', hue='fiber_diam', errorbar='se', palette="RdPu")
# plt.title('Threshold Percent Error by sample and fiber diameter')
plt.legend(title='Fiber Diameter (μm)')
plt.xlabel('Sample')
plt.ylabel('Percent Error')
plt.gca().set_aspect(0.08)
plt.figure()
sns.barplot(data=multimatched, x='nerve_label', y='zdiff_abs', errorbar='se')
plt.title('Activation Z Position Difference by sample')
plt.ylabel('Z Position Difference (mm, absolute value)')
plt.figure()
sns.scatterplot(data=multimatched, x='zdiff_abs', y='pe')
plt.title('Activation Z Position Difference vs. Threshold Percent Error')
plt.ylabel('Percent Error')
plt.xlabel('Z Position Difference (mm, absolute value)')

usedata = multimatched.rename(columns={'threshold': 'threshold2d'})
for comparison in [['pe', 'zdiff'], ['pe', 'zdiff_abs'], ['activation_zpos', 'activation_zpos3d']]:
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
    plt.title(f'Correlation between {comparison[0]} and {comparison[1]}')
#%% percent error of dose response (onset, half, sat, for population and individuals)
#%% Dose-response info all diams
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
plt.figure()
levels = {
    'onset': 0.01,
    'saturation': 1,
}
grouped = threshload.query('sim==3').groupby(['sample', 'fiber_diam', 'sim', 'type', 'nerve_label', 'model', 'nsim'])
analysis = grouped.agg({'threshold': [np.amin, np.median, np.amax]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.rename(columns={'threshold_amin': 'onset', 'threshold_amax': 'saturation'}, inplace=True)
analysis = analysis.reset_index()
# combine onset, saturation, and half into one column with identifier
compiled_data = analysis.melt(
    id_vars=['sample', 'fiber_diam', 'sim', 'type', 'nerve_label', 'model', 'nsim'],
    value_vars=['onset', 'saturation'],
    var_name='level',
    value_name='threshold',
)
# subtract 1 from index if type is 3D
for fiber_diam in pd.unique(compiled_data.fiber_diam):
    ax = sns.boxplot(
        data=compiled_data.query(f'fiber_diam=={int(fiber_diam)}'),
        x='level',
        y='threshold',
        hue='type',
        palette='colorblind',
        color='w',
        boxprops={'facecolor': "none"},
    )
    sns.swarmplot(
        data=compiled_data.query(f'fiber_diam=={int(fiber_diam)}'),
        x='level',
        y='threshold',
        hue='type',
        palette='colorblind',
        dodge=True,
    )
    ax.set_ylim([0, None])
    ax.set_ylabel('Threshold (mA)')
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[2:], labels[2:])
    plt.title(f'fiber_diam={int(fiber_diam)} μm')
    plt.show()
# TODO get this to be paired boxplots again an dline up datapoints
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
g.map_dataframe(sns.boxplot, x='type', y='Threshold (mA)', palette='colorblind', boxprops={'facecolor': "none"})
g.map_dataframe(sns.swarmplot, x='type', y='Threshold (mA)', color='black')
g.map_dataframe(sns.lineplot, x='type', y='Threshold (mA)', units='units', hue='nerve_label', estimator=None, color='k')
plt.subplots_adjust(top=0.9)
plt.legend()
# plt.suptitle('Dose-response curves changes between model types (2D is cathodic leading only)')
# get datamatch for compiled_data
plt.figure()
from src.core.plotter import datamatch_agg

errordr = datamatch_agg(compiled_data.query('type=="2D"'), compiled_data.query('type=="3D"'), 'threshold').drop(
    columns='type'
)
errordr['pe'] = errordr.apply(lambda row: pe(row['threshold3d'], row['threshold']), axis=1)
ax = sns.barplot(data=errordr, hue='fiber_diam', x='level', y='pe', errorbar='se', palette='RdPu')
# plt.title('Percent error between 2D and 3D dose-response (per sample)')
plt.ylabel('Threshold Percent error')
# lineplot
plt.figure()
ax = sns.lineplot(
    data=errordr, hue='sample', x='fiber_diam', style="level", y='pe', errorbar='se', palette='colorblind'
)
plt.ylabel('Threshold Percent error')
# lineplot
plt.figure()
ax = sns.lineplot(data=errordr, x='fiber_diam', hue="level", y='pe', errorbar='se', palette='colorblind')
plt.ylabel('Threshold Percent error')

plt.figure()
# now calculate mean population threshold for 2D and 3D
popdat = errordr.groupby(['fiber_diam', 'level'])['threshold', 'threshold3d'].mean()
popdat.reset_index(inplace=True)
popdat['pe'] = popdat.apply(lambda row: pe(row['threshold3d'], row['threshold']), axis=1)
popdat['level'] = pd.Categorical(popdat['level'], categories=['onset', 'saturation'], ordered=True)
plt.ylim(ax.get_ylim())
sns.barplot(data=popdat, hue='fiber_diam', x='level', y='pe', errorbar='se', palette='RdPu')
# plt.title('Percent error between 2D and 3D dose-response (population mean)')
plt.ylabel('Threshold Percent error')

#%% peri-thk with 2d
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
diam = 3
nsimdata = matched.query(f'fiber_diam=={diam}')
sns.lmplot(data=nsimdata, x='threshold', y='peri_thk', hue='nerve_label', palette='colorblind')
# plot one one line out to max of 3d
plt.ylabel('Perineurium Thickness (micron)')
plt.xlabel('2D Threshold (mA)')
plt.title(f'Fiber Diameter: {diam} μm')
plt.show()

#%% peri-thk with 3d
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
diam = 3
nsimdata = matched.query(f'fiber_diam=={diam}')
sns.lmplot(data=nsimdata, x='threshold3d', y='peri_thk', hue='nerve_label', palette='colorblind')
# plot one one line out to max of 3d
plt.ylabel('Perineurium Thickness (micron)')
plt.xlabel('3D Threshold (mA)')
plt.title(f'Fiber Diameter: {diam} μm')
plt.show()
#%% peri-thk with 2d3d
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
diam = 3
nsimdata = threshload.query(f'fiber_diam=={diam}')
g = sns.lmplot(
    data=nsimdata,
    y='threshold',
    x='peri_thk',
    hue='nerve_label',
    palette='colorblind',
    row='type',
    scatter_kws={'linewidths': 1, 'edgecolor': 'k'},
)
# plot one one line out to max of 3d
g.axes.flat[1].set_xlabel('Threshold (mA)')
for ax in g.axes.flat:
    ax.set_ylabel('Perineurium Thickness (micron)')
# plt.xlabel('3D Threshold (mA)')
# plt.title(f'Fiber Diameter: {diam} μm')
plt.show()
#%% peri-thk with 2d3d
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
diam = 3
nsimdata = threshload.query(f'fiber_diam=={diam}')
sns.lmplot(data=nsimdata, x='threshold', y='peri_thk_act_site', hue='nerve_label', palette='colorblind', col='type')
# plot one one line out to max of 3d
plt.ylabel('Perineurium Thickness (micron)')
# plt.xlabel('3D Threshold (mA)')
# plt.title(f'Fiber Diameter: {diam} μm')
plt.show()
#%% organization compare type
# remove all samples which dont end with a 2
sns.reset_orig()
sns.set(font_scale=1.25, style='whitegrid')
mpl.rcParams['figure.dpi'] = 400
threshdat = threshload[~threshload['sample'].astype(str).str.endswith('0')]
# subtract 1 from sample if type is 3D
threshdat['sample'] = threshdat.apply(lambda x: x['sample'] - 1 if x.type == "3D" else x['sample'], axis=1)
for nerve in pd.unique(threshdat['nerve_label']):
    for n in [0, 5]:
        plt.figure()
        data = threshdat.query(f'nsim=={n} and sim == 3 and nerve_label=="{nerve}"').rename(
            columns={"threshold": "Threshold (mA)", "nerve_label": "Sample", "fiber_diam": "Fiber Diameter (μm)"}
        )
        sns.boxplot(data=data, x='type', y='Threshold (mA)', palette='colorblind', boxprops={'facecolor': "none"})
        sns.stripplot(data=data, jitter=False, linewidth=1, x='type', y='Threshold (mA)', palette='colorblind')
        sns.lineplot(
            data=data,
            x='type',
            y='Threshold (mA)',
            units='master_fiber_index',
            estimator=None,
            color='k',
            # alpha=0.25,
            linewidth=1,
        )
#%% all dose-response fascicle scaled
# # TODO add percent activated to initial processing
# sns.set(font_scale=1.25)
# sns.set_style('whitegrid')
# drdat = threshload.copy()
# drdat['percent_activated'] = 0
# drdat = drdat.rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
# drdat['fiberweight'] = np.nan
# #first need to compute fascicle (inner) weights (percentage of area that inner is in the sample)
# for samp in pd.unique(drdat['sample']):
#     for inner in pd.unique(drdat['inner']):
#         #find the percentage of total peri thk that the inner is
#         total_peri = drdat.query(f'sample=={samp}')['peri_thk'].sum()
#         #check that all inner entries have the same peri thk
#         assert len(pd.unique(drdat.query(f'sample=={samp} and inner=={inner}')['peri_thk'])) == 1
#         percent_peri = drdat.query(f'sample=={samp} and inner=={inner}')['peri_thk'][0]
#         #fiber weight is number of fibers in the inner divided by the total perineurium thickness
