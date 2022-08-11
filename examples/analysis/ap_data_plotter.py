#!/usr/bin/env python3.7

"""The copyrights of this software are owned by Duke University.

Please refer to the LICENSE.txt and README.txt files for licensing
instructions. The source code can be found on the following GitHub
repository: https://github.com/wmglab-duke/ascent
"""

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb

sys.path.append('.')

from src.core.plotter import rename_var
from src.core.query import Query
from src.utils import Object
from src.utils.nd_line import nd_line


#%%
def get_actual_zpos(dat3d, samp3d, model, sim):
    fiberdir = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), '3D_fiberset')
    fibers3d = [x for x in os.listdir(fiberdir) if x.endswith('.dat')]
    for file in fibers3d:
        f_ind = int(os.path.splitext(file)[0])
        coord = np.loadtxt(os.path.join(fiberdir, file), skiprows=1)
        fiberline = nd_line(coord)
        datnew = dat3d[(dat3d["model"] == model) & (dat3d["sim"] == sim) & (dat3d["master_fiber_index"] == f_ind)]
        for index, row in datnew.iterrows():
            actual_zpos = fiberline.interp(row['long_ap_pos'])[2]
            dat3d.loc[index, 'activation_zpos'] = actual_zpos
    return dat3d


cuff_contacts = [21, 29]
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
                'partial_matches': True,
                'include_downstream': True,
                'indices': {'sample': [samp2d], 'model': [0], 'sim': [3]},
            }
        ).run()
        dat2d = q.data()
        dat2d['threed'] = False
        q3 = Query(
            {
                'partial_matches': True,
                'include_downstream': True,
                'indices': {'sample': [253], 'model': [0], 'sim': [3]},
            }
        ).run()
        dat3d = q3.data(source_sample=250)
        dat3d['threed'] = True
        sample_obj = q.get_object(Object.SAMPLE, [250])
        sim_obj = q.get_object(Object.SIMULATION, [250, 0, 3])

        #%%
        dat3z = get_actual_zpos(dat3d, samp3d, model, sim)

        dat2d = dat2d.rename(columns={'long_ap_pos': 'activation_zpos'})

        apdat = pd.concat([dat3z, dat2d])

        #%%
        redict = {"sample": {samp2d: '2D', samp3d: '3D'}}
        datre = rename_var(apdat, redict)
        datre.loc[:, 'activation_zpos'] = datre.loc[:, 'activation_zpos'] / 1000
        #%%
        g = sb.catplot(
            x="sample",
            y='activation_zpos',
            hue='sample',
            col="nsim",
            data=datre,
            kind='strip',
            height=5,
            aspect=0.4,
            linewidth=0,
            order=['2D', '3D'],
            sharey=True,
        )
        axs = g.axes
        axs[0][0].set_ylabel(u'Activation z-position (mm)')
        # for i, s in enumerate([2, 5, 8, 11, 13]):
        #     axs[0][i].set_title(f'fiber diam: {s}μm')
        plt.subplots_adjust(top=0.88)
        plt.suptitle(f"Activation z-position for sample {samp2d}", fontsize=15)
        for ax in axs[0]:
            ln = ax.hlines([cuff_contacts], 0.4, 0.6, color='red')
        axs[0][-1].legend([ln], ['Cuff Contacts'], loc='center right')
        g.savefig(r'out/analysis/{}_zpos.png'.format(samp2d), dpi=400)
