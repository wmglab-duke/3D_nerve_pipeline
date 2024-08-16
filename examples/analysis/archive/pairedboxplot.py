import os

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


# TODO: make loop
os.chdir('../../')
threshdat = pd.read_csv('thresh.csv')
sns.set_style('whitegrid')
# remove all samples which dont end with a 2
threshdat = threshdat[threshdat['sample'].astype(str).str.endswith('2')]
data = split_matched(threshdat).query('nsim == 0 and sim == 3')
# sns.boxplot(data=data, x='type', y='threshold',palette='colorblind')
# sns.stripplot(data=data, x='type', y='threshold',s=10,palette='colorblind')
# sns.lineplot(
#    data=data, x="type", y="threshold", units = 'index',estimator=None,color='k'
# )
# plt.show()
# create facetgrid with sample as columns
g = sns.FacetGrid(data, col="sample", col_wrap=2, sharey=False)
g.map_dataframe(sns.boxplot, x='type', y='threshold', palette='colorblind')
g.map_dataframe(sns.stripplot, x='type', y='threshold', s=10, palette='colorblind')
g.map_dataframe(sns.lineplot, x='type', y='threshold', units='index', estimator=None, color='k')
