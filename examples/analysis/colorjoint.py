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
from scipy.stats import pearsonr

sys.path.append('.')
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
        #%%
        sb.set(font_scale=1.5)
        dat2d = dat2d[dat2d['sample'] == samp2d]
        for nsim in pd.unique(dat2d.nsim):
            plotdata = dat2d.query(f'nsim=={nsim}')
            g = sb.jointplot(data=plotdata, x="threshold3d", y="threshold", hue="inner", palette='colorblind')
            ax = g.ax_joint
            idat = plotdata
            min_thresh = min([min(idat.threshold), min(idat.threshold3d)])
            max_thresh = max([max(idat.threshold), max(idat.threshold3d)])
            limits = (min_thresh, max_thresh)
            ax.set_xlim(limits)
            ax.set_ylim(limits)
            ax.plot(limits, limits, color='red')
            plt.savefig(f'out/analysis/colorjoint{nerve_label}-{samp2d}-{nsim}.png', dpi=500)
