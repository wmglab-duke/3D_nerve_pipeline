#!/usr/bin/env python3.7

"""The copyrights of this software are owned by Duke University.

Please refer to the LICENSE.txt and README.txt files for licensing
instructions. The source code can be found on the following GitHub
repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

import json

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
from scipy.stats import pearsonr

from src.core.plotter import datamatch, rename_var
from src.core.query import Query

#%%
sim = 3
model = 0
with open('examples/analysis/plotconfig.json') as f:
    config = json.load(f)
for sample_data in config['sample_data']:
    samp3d = sample_data['index3d']
    nerve_label = sample_data['name']
    for extrusion_sample in sample_data['exsamples']:
        samp2d = extrusion_sample['index']
        #%%
        q = Query(
            {
                'partial_matches': False,
                'include_downstream': True,
                'indices': {'sample': [samp2d], 'model': [model], 'sim': [sim]},
            }
        ).run()
        dat2d = q.data()

        q = Query(
            {
                'partial_matches': False,
                'include_downstream': True,
                'indices': {'sample': [samp3d], 'model': [model], 'sim': [sim]},
            }
        ).run()

        dat3d = q.data(source_sample=samp2d)

        dat2d = datamatch(dat2d, dat3d, 'threshold')
        #%% Renaming
        redict = {
            # "nsim": {
            #     0: 'fiber diameter: 2\u03BCm',
            #     1: 'fiber diameter: 5\u03BCm',
            #     2: 'fiber diameter: 8\u03BCm',
            #     3: 'fiber diameter: 11\u03BCm',
            #     4: 'fiber diameter: 13\u03BCm',
            # }
        }
        dat2d = dat2d.rename(columns={'threshold': '2D', 'threshold3d': '3D'})
        datre = rename_var(dat2d, redict)
        dat2dnew = datre.drop(columns='3D').rename(columns={'2D': 'threshold'})
        dat2dnew['dataset'] = '2D'
        dat3dnew = datre.drop(columns='2D').rename(columns={'3D': 'threshold'})
        dat3dnew['dataset'] = '3D'
        datfinal = pd.concat([dat2dnew, dat3dnew])

        # datre = dat2d
        #%%
        sb.set(font_scale=1.5)
        plotdata = datfinal[datfinal['sample'] == samp2d]
        g = sb.catplot(
            data=plotdata,
            kind='swarm',
            col='nsim',
            hue='inner',
            y='threshold',
            x='dataset',
            sharey=False,
            palette='colorblind',
        )
        plt.subplots_adjust(top=0.85)
        plt.suptitle(f'Activation thresholds by fascicle (Sample {nerve_label}, 2D slice: {samp2d})')
        axs = g.axes.ravel()
        axs[0].set_ylabel('Activation threshold (mA)')
        plt.subplots_adjust(top=0.85)
        for i, ax in enumerate(g.axes.ravel()):
            corr = {}
            thisdat = dat2d[(dat2d["nsim"] == i) & (dat2d["sample"] == samp2d)]
            corr[samp2d] = round(pearsonr(thisdat['2D'], thisdat['3D'])[0], 3)
            leg = ax.legend(
                labels=["r=" + str(corr[samp2d])], handlelength=0, handletextpad=0, fancybox=True, loc='lower center'
            )
            for item in leg.legendHandles:
                item.set_visible(False)
        plt.savefig(f'out/analysis/colorthresh{nerve_label}-{samp2d}.png', dpi=500)
