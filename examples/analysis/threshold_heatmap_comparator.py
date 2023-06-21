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

from src.core.plotter import heatmaps
from src.core.query import Query

samp2ds = [572, 5729, 5721]
model = 0
simint = 3
samp3ds = [573, 573, 5731]
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
    threshdat.query('nsim in [0,5]', inplace=True)
    threshdat = addpwfd(threshdat, '3')
    # %%
    # calculate activation order
    drdat = threshdat.copy().reset_index()  # change this to repeated?
    drdat['percent_activated'] = 0
    drdat = drdat.rename(columns={'sample': 'samplenum'})
    for i, row in drdat.iterrows():
        thisdat = drdat.query(
            'threed == @row.threed and samplenum == @row.samplenum and fiber_diam == @row.fiber_diam and sim == @row.sim'
        )
        # percent is number of thresholds less than or equal to this threshold divided by total number of thresholds
        drdat.loc[i, 'percent_activated'] = len(thisdat.query('threshold <= @row.threshold')) / len(thisdat)
    drdat = drdat.rename(columns={'samplenum': 'sample'})

    # %%plot heatmaps
    g = sns.FacetGrid(threshdat, row='fiber_diam', col='sample', sharex=False, sharey=False, margin_titles=True)
    g.map(
        heatmaps,
        *threshdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={'s': 25},
        min_max_ticks=True,
    )
    g.set_xlabels('')
    g.set_ylabels('')
    for ax, diam in zip(g.axes[:, 0], [3, 13]):
        ax.set_ylabel(f'D: {diam} um')
    g.set_titles(col_template="col_name", row_template='')

    # plot activation order
    drdat['threshold'] = drdat['percent_activated']
    g = sns.FacetGrid(drdat, row='fiber_diam', col='sample', sharex=False, sharey=False, margin_titles=True)
    g.map(
        heatmaps,
        *drdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={'s': 25},
        cmap=plasmap,
        colorbar=False,
    )
    add_act_colorbar(g.axes)
    g.set_xlabels('')
    g.set_ylabels('')
    for ax, diam in zip(g.axes[:, 0], [3, 13]):
        ax.set_ylabel(f'D: {diam} um')
    for ax, diam in zip(g.axes[:, 0], [3, 13]):
        ax.set_ylabel(f'D: {diam} um')
    g.set_titles(col_template="col_name", row_template='')
