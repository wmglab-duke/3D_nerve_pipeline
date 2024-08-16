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

# samplenames = ['2L', '2R','3R','5R', '6L','6R',  '2Ldef','2Lasc', '3Rdef', '5Rdef', '6Rdef']
# samp2ds = [252, 272,372,572,652,672,2521,2529]
# model = 0
# simint = 10
# samp3ds = [253, 273,373,573,653,673,2531,253]
samplenames = ['5Rdef', '2Ldef', '6Rdef', '3Rdef', '5Rdef']
samp2ds = [5721, 2521, 6721, 3701, 5701]
model = 0
simint = 10
samp3ds = [5731, 2531, 6731, 3731, 5731]
import pickle

import matplotlib.colors as mplcolors
import numpy as np
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
for samp2d, samp3d, samplename in zip(samp2ds, samp3ds, samplenames):
    plotdir = os.path.join(os.getcwd(), fr'plots\activation\{samp2d}-{samp3d}\{simint}')
    os.makedirs(plotdir, exist_ok=True)
    q = Query(
        {
            'partial_matches': True,
            'include_downstream': True,
            'indices': {'sample': [samp2d], 'model': [model], 'sim': [simint]},
        }
    ).run()
    dat2d = q.data(thresh_only=False, source_sim=3)
    dat2d['threed'] = False
    q3 = Query(
        {
            'partial_matches': True,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [simint]},
        }
    ).run()
    dat3d = q3.data(source_sample=samp2d, thresh_only=False, source_sim=3)
    dat3d['threed'] = True
    sample_obj = q.get_object(Object.SAMPLE, [samp2d])
    sim_obj = q.get_object(Object.SIMULATION, [samp2d, model, simint])
    threshdat = pd.concat([dat2d, dat3d])
    threshdat = addpwfd(threshdat, str(simint))
    threshdat.query('fiber_diam in [3,13]', inplace=True)

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
        scatter_kws={'s': 15},
        min_max_ticks=True,
    )
    g.set_xlabels('')
    g.set_ylabels('')
    for ax, diam in zip(g.axes[:, 0], [3, 13]):
        ax.set_ylabel(f'D: {diam} um')
    g.set_titles(col_template="{col_name}", row_template='')
    g.savefig(os.path.join(plotdir, 'threshold.png'), dpi=400)

    # plot activation order
    drdat['threshsave'] = drdat['threshold']
    drdat['threshold'] = drdat['percent_activated']
    g = sns.FacetGrid(drdat, row='fiber_diam', col='sample', sharex=False, sharey=False, margin_titles=True)
    g.map(
        heatmaps,
        *drdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={'s': 15},
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
    g.set_titles(col_template="{col_name}", row_template='')
    g.savefig(os.path.join(plotdir, 'activation.png'), dpi=400)
    # sys.exit()
    # %%
    for plotdiam in [3, 13]:
        # make plot dir
        diamdir = os.path.join(plotdir, f'{plotdiam}', 'activation_error')
        threshdir = os.path.join(plotdir, f'{plotdiam}', 'threshold_error')
        os.makedirs(diamdir, exist_ok=True)
        os.makedirs(threshdir, exist_ok=True)
        # calculate a new column which is error between 2D and 3D activation order
        drdat['activation_error'] = np.nan
        drdat['threshold_error'] = np.nan
        for mfi in drdat['master_fiber_index'].unique():
            for fiber_diam in drdat['fiber_diam'].unique():
                thisdat = drdat.query('master_fiber_index == @mfi and fiber_diam == @fiber_diam')
                thisdat = thisdat.sort_values('threed')
                assert len(thisdat) == 2
                thisdat['activation_error'] = [
                    (np.diff(thisdat['percent_activated'])),
                    (np.diff(thisdat['percent_activated'])),
                ]
                thisdat['threshold_error'] = [(np.diff(thisdat['threshsave'])), (np.diff(thisdat['threshsave']))]
                drdat.loc[thisdat.index, 'activation_error'] = thisdat['activation_error']
                drdat.loc[thisdat.index, 'threshold_error'] = thisdat['threshold_error']
        fiberpath = os.path.join(os.getcwd(), fr'samples\{samp3d}\models\0\sims\3\3D_fiberset')
        slidespath = os.path.join(os.getcwd(), fr'input\slides\{samplename}slides.obj')

        # load pickled slidelist
        with open(slidespath, 'rb') as f:
            slidelist = pickle.load(f)

        fibers = {}
        # load each fiber file and append to list
        for file in os.listdir(fiberpath):
            # print(file)
            # if file == '1.dat':# use this to watch a specific fiber
            if file.endswith('.dat'):
                fibers[int(file.replace('.dat', ''))] = np.loadtxt(os.path.join(fiberpath, file), skiprows=1)

        # %%
        # loop through each slide, and calculate z position
        for i, slide in enumerate(slidelist):
            if i % 5 != 0:  # every 100 microns
                continue
            # slide.scale(0.5)
            zpos = i * 20  # 20 microns per slice
            # contact centers are 29050 and 21050
            # if zpos < 20000 or zpos > 30000:
            #     continue  # use this to check only a specific range
            if zpos < 27000 or zpos > 31000:
                continue  # use this to check only a specific range
            # get list of x,y coordinates for each fiber at this z position
            xs = []
            ys = []
            xsact = []
            ysact = []
            colors = []
            colorsthresh = []
            for mfi, fiber in fibers.items():
                # find the index of the closest z value
                idx = (np.abs(fiber[:, 2] - zpos)).argmin()
                # get the x,y coordinates at that index
                x = fiber[idx, 0]
                y = fiber[idx, 1]
                # append to slide list
                xs.append(x)
                ys.append(-y)
                zposact = drdat.query('master_fiber_index==@mfi and fiber_diam==@plotdiam and sample==@samp3d')[
                    'activation_zpos'
                ].values[0]
                if np.abs(zpos - zposact) <= 50:
                    xsact.append(x)
                    ysact.append(-y)
                # color according to threshold value, using plasmap
                # color = plasmap(drdat.query('master_fiber_index==@mfi and fiber_diam==@plotdiam and sample==@samp3d')['threshold'].values[0])
                # colors.append(color)
                # now color according to error, but use a colormap with 0 as white
                # # newcolormap = colormaps['binary']
                # newcolormap = colormaps['PiYG']
                # color = newcolormap(drdat.query('master_fiber_index==@mfi and fiber_diam==@plotdiam and sample==@samp3d')['error'].values[0])
                # colors.append(color)
                colors.append(
                    drdat.query('master_fiber_index==@mfi and fiber_diam==@plotdiam and sample==@samp3d')[
                        'activation_error'
                    ].values[0]
                )
                colorsthresh.append(
                    drdat.query('master_fiber_index==@mfi and fiber_diam==@plotdiam and sample==@samp3d')[
                        'threshold_error'
                    ].values[0]
                )
            # plot the slide and all fiber points
            plt.figure()
            plt.scatter(xs, ys, s=3, c=colors, cmap=colormaps['PiYG'], vmin=-0.5, vmax=0.5)
            plt.colorbar()
            plt.scatter(xsact, ysact, c='k', marker='x')
            # label the top of colorbar "activates first in 2D" and bottom "activates first in 3D"
            plt.text(1.25, 1.03, 'activates first in 2D', ha='center', va='bottom', transform=plt.gca().transAxes)
            plt.text(1.25, -0.03, 'activates first in 3D', ha='center', va='top', transform=plt.gca().transAxes)
            plt.title(f'Slide {i}-zpos{zpos}')
            # add colorbar with new colormap
            slide.plot(final=False)
            # plt.show()
            plt.savefig(os.path.join(diamdir, f'slide{i}_zpos{zpos}.png'))

            # now figure for threshold
            plt.figure()
            maxval = np.amax([np.abs(x) for x in colorsthresh])
            plt.scatter(xs, ys, s=3, c=colorsthresh, cmap=colormaps['PiYG'], vmin=-maxval, vmax=maxval)
            plt.colorbar()
            plt.scatter(xsact, ysact, c='k', marker='x')
            # label the top of colorbar "higher threshold in 2D" and bottom "higher threshold in 3D"
            plt.text(1.25, 1.03, 'lower threshold in 2D', ha='center', va='bottom', transform=plt.gca().transAxes)
            plt.text(1.25, -0.03, 'lower threshold in 3D', ha='center', va='top', transform=plt.gca().transAxes)
            plt.title(f'Slide {i}-zpos{zpos}')
            # add colorbar with new colormap
            slide.plot(final=False)
            # plt.show()
            plt.savefig(os.path.join(threshdir, f'slide{i}_zpos{zpos}_thresh.png'))
