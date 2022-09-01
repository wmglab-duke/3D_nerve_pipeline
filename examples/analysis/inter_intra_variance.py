#!/usr/bin/env python3
import json
import os
import sys

import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.append('.')

from src.core.plotter import get_datamatch
from src.core.query import Query
from src.utils import Object

model = 0
with open('examples/analysis/plotconfig.json') as f:
    config = json.load(f)
for simdex in config['sim_data'].keys():
    simint = int(simdex)
    for sample_data in config['sample_data']:
        samp3d = sample_data['index3d']
        nerve_label = sample_data['name']
        samples2d = [x['index'] for x in sample_data['exsamples']]
        #%%
        threshdat = get_datamatch(samples2d, samp3d, model, simint, nerve_label)
        test = threshdat.groupby(['inner', 'sample', 'nsim'], as_index=False)
        analysis = test.agg({'threshold': [np.var, np.mean], 'threshold3d': [np.var, np.mean]})
        analysis.columns = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
        analysis3d = analysis.copy()
        analysis.drop('threshold3d_mean', axis=1, inplace=True)
        analysis.drop('threshold3d_var', axis=1, inplace=True)
        analysis3d.drop('threshold_mean', axis=1, inplace=True)
        analysis3d.drop('threshold_var', axis=1, inplace=True)
        analysis3d.rename(
            columns={'threshold3d_mean': 'threshold_mean', 'threshold3d_var': 'threshold_var'}, inplace=True
        )
        analysis['threed'] = False
        analysis3d['threed'] = True
        data = pd.concat([analysis, analysis3d])
        data = data.query('nsim==0')
        # variance within inners (lower in 2d)
        sns.catplot(
            data=data,
            kind='strip',
            x='threed',
            hue='inner',
            y='threshold_var',
            col='sample',
            sharex=False,
            sharey=False,
            dodge=True,
        )
        sns.catplot(
            data=data,
            kind='strip',
            x='threed',
            hue='inner',
            y='threshold_mean',
            col='sample',
            sharex=False,
            sharey=False,
            dodge=True,
        )
        sns.boxplot(data=data, hue='threed', y='threshold_var', x='sample')
        plt.title('within-fascicle threshold variance')
        plt.figure()
        sns.boxplot(data=data, hue='threed', y='threshold_mean', x='sample')
        plt.title('within-fascicle threshold mean')
