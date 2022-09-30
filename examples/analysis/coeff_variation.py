import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import variation

os.chdir('../../')


threshdat = pd.read_csv('thresh.csv')
# melt threshdat threshold and threshold3d to one column
threshdat.rename(columns={'threshold3d': '3D', 'threshold': '2D'}, inplace=True)
threshdat = threshdat.melt(
    id_vars=['sample', 'nsim', 'sim', 'master_fiber_index', 'inner'],
    value_vars=['2D', '3D'],
    var_name='type',
    value_name='threshold',
)
sns.catplot(
    data=threshdat.query('nsim in [0,5]'),
    x='sample',
    y='threshold',
    hue='type',
    kind='box',
    palette='colorblind',
    col='nsim',
    sharey=False,
)
plt.suptitle('2D vs 3D Thresholds for all samples')

grouped = threshdat.groupby(['sample', 'nsim', 'sim', 'type'])
# calculate coefficient of variation for thresholds on grouped data
analysis = grouped.agg({'threshold': [variation]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis = analysis.reset_index()
sns.catplot(data=analysis, x='type', y='threshold_variation', hue='nsim', kind='box', palette='colorblind')
plt.ylabel('Threshold Coefficient of Variation')
plt.title('Coefficient of Variation of Thresholds for 2D and 3D')


vardat = threshdat.query('nsim in [0,5] and sim==3')
grouped = vardat.groupby(['sample', 'nsim', 'sim', 'type', 'inner'])
analysis = grouped.agg({'threshold': [np.var, np.mean]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.reset_index(inplace=True)
sns.catplot(
    data=analysis, sharey=False, col='nsim', x='sample', y='threshold_var', hue='type', kind='bar', palette='colorblind'
)
plt.suptitle('Variance values of thresholds within fascicles')
sns.catplot(
    data=analysis, sharey=False, x='type', col='nsim', y='threshold_var', kind='bar', palette='colorblind', aspect=0.5
)
plt.subplots_adjust(top=0.9)
plt.suptitle('Variance values of thresholds within fascicles')

# now do variance between fascicle mean thresholds
grouped = analysis.groupby(['sample', 'nsim', 'type'])
analysis = grouped.agg({'threshold_mean': [np.var]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.reset_index(inplace=True)
sns.catplot(
    data=analysis,
    sharey=False,
    col='nsim',
    x='sample',
    y='threshold_mean_var',
    hue='type',
    kind='bar',
    palette='colorblind',
)
plt.suptitle('Variance values of fascicle thresholds')
sns.catplot(
    data=analysis,
    sharey=False,
    x='type',
    col='nsim',
    y='threshold_mean_var',
    kind='bar',
    palette='colorblind',
    aspect=0.5,
)
plt.subplots_adjust(top=0.9)
plt.suptitle('Variance values of fascicle thresholds')


threshdat['z_score'] = 0
threshdat['z_score_error'] = 0
threshdat = threshdat.rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
for i, row in threshdat.iterrows():
    thisdat = threshdat.query(
        'modeltype == @row.modeltype and samplenum == @row.samplenum and nsim == @row.nsim and sim == @row.sim'
    )
    # percent is number of thresholds less than or equal to this threshold divided by total number of thresholds
    threshdat.loc[i, 'z_score'] = (row.threshold - np.mean(thisdat.threshold)) / np.std(thisdat.threshold)
for i, row in threshdat.iterrows():
    # subtract this row z score from the z score of nsim 0
    threshdat.loc[i, 'z_score_error'] = float(
        row.z_score
        - threshdat.query(
            'master_fiber_index == @row.master_fiber_index and modeltype == @row.modeltype and samplenum == @row.samplenum and nsim == 0 and sim == @row.sim'
        ).z_score
    )
sns.catplot(data=threshdat, x='modeltype', y='z_score_error', hue='nsim', kind='box', palette='colorblind')
plt.title('Z Score Residual for 2D and 3D')
plt.ylabel('Z Score Residual from nsim=0')
