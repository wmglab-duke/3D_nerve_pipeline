#!/usr/bin/env python3.7

"""The copyrights of this software are owned by Duke University.

Please refer to the LICENSE.txt and README.txt files for licensing
instructions. The source code can be found on the following GitHub
repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

import json
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('agg')

sys.path.append('.')

from src.core.plotter import corrcalc, get_datamatch

model = 0

with open('examples/analysis/plotconfig.json') as f:
    config = json.load(f)


def run_comparison(comparison):
    allcorr = []
    plt.figure()
    simint = int(simdex)
    os.makedirs(f'out/analysis/{simdex}', exist_ok=True)
    for sample_data in config['sample_data']:
        samp3d = sample_data['index3d']
        nerve_label = sample_data['name']
        samples2d = [x['index'] for x in sample_data['exsamples']]
        # PLOT CORRELATION
        inputdata = get_datamatch(
            samples2d,
            samp3d,
            model,
            simint,
            nerve_label,
            tortuosity='tortuosity3d' in comparison,
            source_sim=config['source_sim'],
        )
        if 'tortuosity3d' in comparison:
            inputdata = inputdata.query(f'sample=={pd.unique(inputdata["sample"])[0]}')
        corrdata = corrcalc(inputdata, comparison)
        for i, d in enumerate(corrdata):
            for k, v in d.items():
                allcorr.append(
                    {
                        "sample": nerve_label,
                        "correlation": v,
                        "contact": "rostral" if str(k).endswith('0') else "caudal",
                        "nsim": i,
                    }
                )
    # plot correlation aggreggate
    plt.figure()
    finalcorr = pd.DataFrame(allcorr)
    sns.swarmplot(data=finalcorr, x='nsim', y='correlation', hue='sample', s=10, palette='colorblind')
    sns.lineplot(
        data=finalcorr, x='nsim', y='correlation', style='contact', hue='sample', legend=False, palette='colorblind'
    )
    plt.title(f'{comparison[0]}-{comparison[1]}')
    plt.gcf().savefig(f'out/analysis/{comparison[0]}-{comparison[1]}-sim{simint}')


for simdex in config['sim_data'].keys():
    plt.close('all')
    run_comparison(['threshold', 'threshold3d'])
    run_comparison(['threshold', 'peri_thk'])
    run_comparison(['threshold3d', 'peri_thk'])
    run_comparison(['threshold3d', 'tortuosity3d'])
