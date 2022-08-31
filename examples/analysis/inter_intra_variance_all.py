#!/usr/bin/env python3
import json
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.append('.')

from src.core.plotter import get_datamatch

datas = []
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
        analysis['model'] = '2D'
        analysis3d['model'] = '3D'
        data = pd.concat([analysis, analysis3d])
        data['nerve_label'] = nerve_label
        datas.append(data)
    alldat = pd.concat(datas)
    sns.catplot(
        data=alldat,
        kind='box',
        col='nsim',
        row='nerve_label',
        hue='model',
        y='threshold_var',
        x='sample',
        sharey=False,
        margin_titles=True,
    )
    plt.suptitle('within-fascicle threshold variance')
    plt.gcf().savefig('fascvar.png', dpi=400, bbox_inches='tight')
    plt.figure()
    sns.catplot(
        data=alldat,
        kind='box',
        col='nsim',
        row='nerve_label',
        hue='model',
        y='threshold_mean',
        x='sample',
        sharey=False,
        margin_titles=True,
    )
    plt.suptitle('within-fascicle threshold mean')
    plt.gcf().savefig('fascmean.png', dpi=400, bbox_inches='tight')
