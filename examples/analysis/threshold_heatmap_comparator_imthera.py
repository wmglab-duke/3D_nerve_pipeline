#!/usr/bin/env python3.7

"""Generate a heatmap of activation thresholds.

The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent.

Note: if more than one heatmap is desired, you must use a Seaborn FacetGrid.
RUN THIS FROM REPOSITORY ROOT
"""
import os

os.chdir('../..')
import json

import matplotlib.pyplot as plt
import numpy as np

from src.core.plotter import heatmaps
from src.core.query import Query

samp2ds = [5715]
model = 0
simint = 17
samp3ds = [5735]
import matplotlib.colors as mplcolors
import pandas as pd
import seaborn as sns
from matplotlib import colormaps

from src.utils import Object


def add_act_colorbar(ax):
    cmap = colormaps['plasma']
    cmap.set_bad(color='w')
    cmap = cmap.reversed()
    mappable = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=mplcolors.Normalize(vmin=0, vmax=1),
    )
    cb = plt.colorbar(mappable=mappable, ax=ax, ticks=[0, 1])
    cb.ax.set_ylabel('Percent Activated')


def add_thresh_colorbar(ax, mint, maxt, **kwargs):
    cmap = colormaps['viridis']
    cmap.set_bad(color='w')
    cmap = cmap.reversed()
    mappable = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=mplcolors.Normalize(vmin=mint, vmax=maxt),
    )
    cb = plt.colorbar(mappable=mappable, ax=ax, ticks=[mint, maxt], **kwargs)
    cb.ax.set_ylabel('Threshold (mA)')


def addpwfd(data, sim):
    with open('examples/analysis/plotconfig.json') as f:
        config = json.load(f)
    nsim_key = config['sim_data'][sim]['nsim_key']
    for nsim in nsim_key:
        pulse_width = nsim_key[nsim]['pulse_width']
        fiber_diam = nsim_key[nsim]['fiber_diam']
        nsim = int(nsim)
        data.loc[data['nsim'] == nsim, 'pulse_width'] = pulse_width
        data.loc[data['nsim'] == nsim, 'fiber_diam'] = fiber_diam
    return data


plasmap = colormaps['plasma']
plasmap.set_bad(color='w')
plasmap = plasmap.reversed()
# %%
for samp2d, samp3d in zip(samp2ds, samp3ds):
    q = Query(
        {
            'partial_matches': True,
            'include_downstream': True,
            'indices': {'sample': [samp2d], 'model': [model], 'sim': [simint]},
        }
    ).run()
    dat2d = q.data(thresh_only=True)
    dat2d['threed'] = False
    q3 = Query(
        {
            'partial_matches': True,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [simint]},
        }
    ).run()
    dat3d = q3.data(source_sample=samp2d, thresh_only=True)
    dat3d['threed'] = True
    sample_obj = q.get_object(Object.SAMPLE, [samp2d])
    sim_obj = q.get_object(Object.SIMULATION, [samp2d, model, simint])
    threshdat = pd.concat([dat2d, dat3d])
    threshdat.query('fiberset_index==0', inplace=True)
    # %%
    # # calculate activation order
    # drdat = threshdat.copy().reset_index()  # change this to repeated?
    # drdat['percent_activated'] = 0
    # drdat = drdat.rename(columns={'sample': 'samplenum'})
    # for i, row in drdat.iterrows():
    #     thisdat = drdat.query(
    #         'threed == @row.threed and samplenum == @row.samplenum and fiber_diam == @row.fiber_diam and sim == @row.sim'
    #     )
    #     # percent is number of thresholds less than or equal to this threshold divided by total number of thresholds
    #     drdat.loc[i, 'percent_activated'] = len(thisdat.query('threshold <= @row.threshold')) / len(thisdat)
    # drdat = drdat.rename(columns={'samplenum': 'sample'})

    # %%plot heatmaps
    g = sns.FacetGrid(threshdat, row='sample', col='active_src_index', sharex=False, sharey=False, margin_titles=True)
    g.map(
        heatmaps,
        *threshdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={'s': 25},
        min_max_ticks=True,
        min_thresh=threshdat.threshold.min(),
        max_thresh=threshdat.threshold.max(),
        colorbar=False,
    )
    minsave = threshdat.threshold.min()
    maxsave = threshdat.threshold.max()
    g.set_xlabels('')
    g.set_ylabels('')
    g.set_titles(col_template="{col_name}", row_template='')
    for ax, name in zip(g.axes[:, 0], ['2DEM', '3DM']):
        # TODO get title and infer name from whether there is a 3
        ax.set_ylabel(name)
    add_thresh_colorbar(g.axes, threshdat.threshold.min(), threshdat.threshold.max(), panchor=(-0.1, 1))
    plt.subplots_adjust(hspace=0, wspace=-0.5)

    for ind, ax in enumerate(g.axes.flat):
        # add a circular band around the nerve with radius 1500 um from 0 to 270 degrees
        ax.set_xlim([val * 1.3 for val in ax.get_xlim()])
        ax.set_ylim([val * 1.3 for val in ax.get_ylim()])

        ax.set_autoscale_on(False)

        theta = np.linspace(0, np.radians(325), 100)
        r = 1700
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        # add 6 evenly spaced stars from 30 to 300 degrees, add numbers 0-5 to the stars
        theta = np.linspace(np.radians(50), np.radians(325 - 50), 6)
        r = 1550
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        i = ind % 6
        ax.plot(x[i], y[i], 'ro', markersize=15, color='r', label='contacts')
        # ax.text(x[i], y[i], str(i), fontsize=10, color='k', ha='center', va='center')
