#!/usr/bin/env python3.7

"""The copyrights of this software are owned by Duke University.

Please refer to the LICENSE.txt and README.txt files for licensing
instructions. The source code can be found on the following GitHub
repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

import json
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

sys.path.append('.')
from src.core.query import Query

#%%
sim = 3
model = 0
with open('examples/analysis/plotconfig.json') as f:
    config = json.load(f)

for sample_data in config['sample_data']:
    samp3d = sample_data['index3d']
    nerve_label = sample_data['name']
    samples2d = [x['index'] for x in sample_data['exsamples']]
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': samples2d, 'model': [model], 'sim': [sim]},
        }
    ).run()

    # builds heatmaps
    # q.barcharts_compare_models(logscale=False,
    #                            model_labels=['Model 0: Veltink Epineurium, \n              Veltink Perineurium',
    #                                          'Model 1: Veltink Epineurium, \n              Goodall Perineurium',
    #                                          'Model 2: Goodall Epineurium, \n              Veltink Perineurium',
    #                                          'Model 3: Goodall Epineurium, \n              Goodall Perineurium']
    #                            )
    dat2d = q.data()

    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [sim]},
        }
    ).run()

    dat3d = q.data(source_sample=samples2d[0])

    data = pd.concat([dat2d, dat3d])
    data.reset_index(inplace=True)

    for i in range(len(pd.unique(data['nsim']))):
        # ax = axs[i]
        sb.color_palette("tab10")
        # plt.figure()
        plotdata = data[data.nsim == i]
        plotdata = plotdata[plotdata['sample'] != 670]
        sb.ecdfplot(data=plotdata, x='threshold', hue='sample', palette='colorblind')
        plt.xscale('log')
        plt.ylabel('Proportion of Fibers Activated')
        plt.xlabel('Activation Threshold (mA, log scale)')
        plt.title(f'Threshold eCDF for {nerve_label}')
    # plt.text(0.05, -0.25, 'Note: fiber diameters (\u03bcm) from left to right: [13, 11, 8, 5, 2]', fontstyle='italic')
    plt.savefig(r'out/analysis/{}_ecdf.png'.format(nerve_label), dpi=400, bbox_inches='tight')
