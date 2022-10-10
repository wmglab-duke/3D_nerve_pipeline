# -*- coding: utf-8 -*-
import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import pearsonr


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


# TODO: make loop
os.chdir('../../')
threshdat = pd.read_csv('thresh.csv')
sns.set_style('whitegrid')
# remove all samples which dont end with a 2
threshdat = threshdat[threshdat['sample'].astype(str).str.endswith('2')]
data = split_matched(threshdat).query('nsim in [0,5] and sim == 3')
# apply minmax normalization within sample and nsim
data['threshold'] = data.groupby(['sample', 'nsim', 'type'])['threshold'].transform(
    lambda x: (x - x.min()) / (x.max() - x.min())
)
# set nsim as categorical
data['nsim'] = data['nsim'].astype('category')
# nsim 0 is three and nsim 5 is thirteen
data['nsim'].cat.rename_categories({0: '3', 5: '13'}, inplace=True)
# g = sns.FacetGrid(data, row="type", col="sample", sharey=False)
# # g.map_dataframe(sns.boxplot, x='nsim', y='threshold', palette='colorblind')
# # g.map_dataframe(sns.stripplot, x='nsim', y='threshold', s=10, palette='colorblind')
# g.map_dataframe(sns.lineplot, x='nsim', y='threshold', units='index', estimator=None, color='k')
g = sns.FacetGrid(data, col="sample", row='type', sharey=False, margin_titles=True)
g.map_dataframe(sns.boxplot, x='nsim', y='threshold', palette='colorblind', boxprops={'facecolor': "none"})
g.map_dataframe(sns.stripplot, linewidth=1, x='nsim', y='threshold', palette='colorblind')
g.map_dataframe(sns.lineplot, x='nsim', y='threshold', units='master_fiber_index', estimator=None, color='k')
plt.subplots_adjust(top=0.9)
plt.suptitle('min-max normalized thresholds compared between 3 um and 13 um thresholds')
plt.savefig('matchsim.png', dpi=400)
