#!/usr/bin/env python3
import json
import sys

import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src.core.plotter import heatmaps
from src.core.query import Query
from src.utils import Object

sys.path.append('.')

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
                'indices': {'sample': [samp2d], 'model': [model], 'sim': [sim]},
            }
        ).run()
        dat2d = q.data()
        dat2d['threed'] = False
        q3 = Query(
            {
                'partial_matches': True,
                'include_downstream': True,
                'indices': {'sample': [samp3d], 'model': [model], 'sim': [sim]},
            }
        ).run()
        dat3d = q3.data(source_sample=samp2d)
        dat3d['threed'] = True
        sample_obj = q.get_object(Object.SAMPLE, [samp2d])
        sim_obj = q.get_object(Object.SIMULATION, [samp2d, model, sim])
        threshdat = pd.concat([dat2d, dat3d])
        #%%
        g = sns.FacetGrid(threshdat, row='nsim', col='sample', sharex=False, sharey=False)
        g.map(heatmaps, *threshdat.columns, sample_object=sample_obj, sim_object=sim_obj, scatter_kws={'s': 25})

        # Title and clear axis labels
        plt.subplots_adjust(top=0.95)
        plt.suptitle('Grid of activation threshold heatmaps')
        for ax in g.axes.ravel():
            ax.set_xlabel('')
            ax.set_ylabel('')
        plt.savefig(f'out/analysis/heatmap3D-{samp2d}-{sim}', dpi=400, bbox_inches='tight')
        #%% Process seaborn data

        def add_colorbar(ax):
            cmap = plt.cm.get_cmap('viridis')
            cmap.set_bad(color='w')
            cmap = cmap.reversed()
            mappable = plt.cm.ScalarMappable(
                cmap=cmap,
                norm=mplcolors.Normalize(vmin=min_prenorm, vmax=max_prenorm),
            )
            cb_label = r'mA'
            cb = plt.colorbar(mappable=mappable, ax=ax)
            cb.ax.set_title(cb_label)

        for n in pd.unique(threshdat.nsim):
            this = threshdat.query(f"nsim=={n}")
            min_prenorm, max_prenorm = min(this.threshold), max(this.threshold)
            this.threshold = this.threshold / max(this.threshold)
            g = sns.FacetGrid(this, col='sample', sharex=False, sharey=False)
            g.map(
                heatmaps,
                *threshdat.columns,
                sample_object=sample_obj,
                sim_object=sim_obj,
                colorbar=False,
                min_thresh=min(this.threshold),
                max_thresh=max(this.threshold),
                scatter_kws={'s': 25},
            )

            # Title and clear axis labels
            plt.subplots_adjust(top=0.8)
            plt.suptitle(f'Paired heatmaps for nsim {n}')
            for ax in g.axes.ravel():
                ax.set_xlabel('')
                ax.set_ylabel('')
            add_colorbar(g.axes)
            plt.savefig(f'out/analysis/heatmap3d-{samp2d}-{sim}-{n}', dpi=400, bbox_inches='tight')
