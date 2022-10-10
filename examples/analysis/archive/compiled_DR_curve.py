# -*- coding: utf-8 -*-
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# TODO loop over sim and nsim


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
# for each sample, nsim, and sim, calculate the percent fibers activated for each threshold
threshdat['percent_activated'] = 0
infodat = split_matched(threshdat)
infodat = infodat.rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
selnsim = 0
infodat = infodat.query(f'sim==3 and nsim=={selnsim}')
for i, row in infodat.iterrows():
    thisdat = infodat.query(
        'modeltype == @row.modeltype and samplenum == @row.samplenum and nsim == @row.nsim and sim == @row.sim'
    )
    # percent is number of thresholds less than or equal to this threshold divided by total number of thresholds
    infodat.loc[i, 'percent_activated'] = len(thisdat.query('threshold <= @row.threshold')) / len(thisdat)
    # could shade in min and max to show response range
sns.scatterplot(data=infodat, x='percent_activated', y='threshold', hue='modeltype', palette='colorblind')
plt.xlabel('Percent of fibers activated')
plt.ylabel('Threshold (mA)')
plt.title(f'Dose-Response for all samples nsim={selnsim}')
