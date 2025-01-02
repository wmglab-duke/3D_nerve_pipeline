#!/usr/bin/env python3
import json
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.append('.')

from src.core.plotter import get_datamatch

alldat = pd.read_csv('thresh.csv')

datas = []
model = 0
with open('examples/analysis/plotconfig.json') as f:
    config = json.load(f)
for simdex in ['3']:
    simint = int(simdex)
    for sample_data in config['sample_data']:
        samp3d = sample_data['index3d']
        nerve_label = sample_data['name']
        samples2d = [x['index'] for x in sample_data['exsamples']]
        # remove all integer samples ending in 0 from samples2d
        samples2d = [x for x in samples2d if not str(x).endswith('0')]
        # %%
        # threshdat query from alldat, sample should be in samples2d and sim should be simdex
        threshdat = alldat[(alldat['sample'].isin(samples2d)) & (alldat['sim'] == simint)]
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
    maxnsim = alldat['nsim'].max()
    alldat = alldat.query(f'nsim in [0,{maxnsim}]')
    sns.catplot(
        data=alldat,
        kind='bar',
        col='nsim',
        # col_wrap=3,
        hue='model',
        y='threshold_var',
        x='sample',
        sharey=False,
        margin_titles=True,
    )
    plt.show()
    plt.suptitle(f'within-fascicle threshold variance sim {simint}')
    plt.gcf().savefig('fascvar.png', dpi=400, bbox_inches='tight')
    plt.figure()
    sns.catplot(
        data=alldat,
        kind='box',
        col='nsim',
        # col_wrap=3,
        hue='model',
        y='threshold_mean',
        x='sample',
        sharey=False,
        margin_titles=True,
    )
    plt.suptitle(f'within-fascicle threshold mean sim {simint}'),
    plt.gcf().savefig('fascmean.png', dpi=400, bbox_inches='tight')
    grouped = alldat.groupby(['sample', 'nsim', 'model'], as_index=False)
    all2 = grouped.agg({'threshold_mean': [np.var]})
    all2.columns = ["_".join(col_name).rstrip('_') for col_name in all2.columns]
    plt.figure()
    sns.catplot(
        data=all2,
        kind='swarm',
        col='nsim',
        hue='sample',
        y='threshold_mean_var',
        x='model',
        sharey=False,
        margin_titles=True,
    )
    plt.suptitle(f'between-fascicle threshold variance sim {simint}')
    plt.gcf().savefig('compilevar.png', dpi=400, bbox_inches='tight')
