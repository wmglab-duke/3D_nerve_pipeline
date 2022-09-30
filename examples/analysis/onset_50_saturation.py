# -*- coding: utf-8 -*-
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def split_matched(data):  # todo replace with pd.melt
    """Split matched data into 2d and 3d dataframes."""
    data2d = data.drop(columns=['threshold3d'])
    data2d['type'] = '2D'
    data3d = data.drop(columns=['threshold'])
    data3d['type'] = '3D'
    data3d.reset_index(inplace=True, drop=False)
    data2d.reset_index(inplace=True, drop=False)
    data3d.rename(columns={'threshold3d': 'threshold'}, inplace=True)
    return pd.concat([data2d, data3d], axis=0, ignore_index=True)


os.chdir('../../')
threshdat = pd.read_csv('thresh.csv')
sns.set_style('whitegrid')
levels = {
    'onset': 0.01,
    'half': 0.5,
    'saturation': 1,
}
# TODO: make loop and save figs and stop speccing nsim
grouped = split_matched(threshdat).query('sim==3 and nsim in [0,5]').groupby(['sample', 'nsim', 'sim', 'type'])
analysis = grouped.agg({'threshold': [np.amin, np.median, np.amax]})
analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
analysis.rename(
    columns={'threshold_amin': 'onset', 'threshold_median': 'half', 'threshold_amax': 'saturation'}, inplace=True
)
analysis = analysis.reset_index()
# combine onset, saturation, and half into one column with identifier
compiled_data = analysis.melt(
    id_vars=['sample', 'nsim', 'sim', 'type'],
    value_vars=['onset', 'half', 'saturation'],
    var_name='level',
    value_name='threshold',
)
# subtract 1 from index if type is 3D
for nsim in pd.unique(compiled_data.nsim):
    ax = sns.boxplot(
        data=compiled_data.query(f'nsim=={int(nsim)}'),
        x='level',
        y='threshold',
        hue='type',
        palette='colorblind',
        color='w',
        boxprops={'facecolor': "none"},
    )
    sns.swarmplot(
        data=compiled_data.query(f'nsim=={int(nsim)}'),
        x='level',
        y='threshold',
        hue='type',
        palette='colorblind',
        dodge=True,
    )
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[2:], labels[2:])
    plt.title(f'nsim={int(nsim)}')
    plt.show()
# TODO get this to be paired boxplots
# set up facetgrid with nsim as row and level as columns
compiled_data.reset_index(inplace=True)
# subtract 1 from sample where type is 3D
compiled_data['index'] = compiled_data.apply(lambda x: x['index'] - 1 if x['type'] == '3D' else x['index'], axis=1)
g = sns.FacetGrid(compiled_data, col="level", row="nsim", sharey=False, margin_titles=True)
g.map_dataframe(sns.boxplot, x='type', y='threshold', palette='colorblind', boxprops={'facecolor': "none"})
g.map_dataframe(sns.swarmplot, x='type', y='threshold', palette='colorblind', dodge=True)
g.map_dataframe(sns.lineplot, x='type', y='threshold', units='index', estimator=None, color='k')
