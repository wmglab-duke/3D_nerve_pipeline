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

import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colormaps
from src.core.plotter import heatmaps
from src.core.query import Query
from src.utils import Object

samp2ds = [571, 5719, 5711]
model = 0
simint = 3
samp3ds = [573, 573, 5731]
pairnames = ["Undeformed", "2D-3D", "3D-3D"]
slidenames = ['5RDS5', '5RdefDS5', '5RdefDS5']


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


def add_thresh_colorbar(ax, mint, maxt):
    cmap = colormaps['viridis']
    cmap.set_bad(color='w')
    cmap = cmap.reversed()
    mappable = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=mplcolors.Normalize(vmin=mint, vmax=maxt),
    )
    cb = plt.colorbar(mappable=mappable, ax=ax, ticks=[mint, maxt])
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
for samp2d, samp3d, pairname, slidename in zip(samp2ds, samp3ds, pairnames, slidenames):
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
    threshdat.threshold /= threshdat.threshold.max()
    threshdat.query('nsim in [0,5]', inplace=True)
    threshdat = addpwfd(threshdat, '3')
    threshdat.query('nsim==0', inplace=True)
    # %%plot heatmaps
    titles = [f'2DEM\n{pairname}', '3DM']
    g = sns.FacetGrid(threshdat, row='fiber_diam', col='sample',
                      sharex=False, sharey=False, margin_titles=True)
    g.map(
        heatmaps,
        *threshdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={'s': 10},
        min_max_ticks=True,
        min_thresh=threshdat.threshold.min(),
        max_thresh=threshdat.threshold.max(),
        colorbar=False
    )
    minsave = threshdat.threshold.min()
    maxsave = threshdat.threshold.max()
    g.set_xlabels('')
    g.set_ylabels('')
    g.set_titles(col_template="{col_name}", row_template='')
    for ax in g.axes[0, :]:
        thename = titles[0] if ax.get_title()[2] != '3' else titles[1]
        ax.set_title(thename)
    add_thresh_colorbar(
        g.axes,
        threshdat.threshold.min(),
        threshdat.threshold.max(),
    )
    # load contact coordinates
    dire = rf'D:\threed_ascent\input\contact_coords\{slidename}'
    for i in range(2):
        contact_coords = np.loadtxt(
            dire + r'\pcs' + str(i + 1) + '.txt', skiprows=8)[:, :2]
        for ax in g.axes.ravel():
            ax.scatter(contact_coords[:, 0], contact_coords[:, 1],
                       s=1, color='k', label='contacts' if i == 0 else '_')
calculate earth mover distance
    import cv2
    from cv2 import EMD
    fiberpoints = np.array([x['fiber'][0][:2]
                           for x in sim_obj.fibersets[0].fibers])
    fiberpoints[:, 1] += -np.amin(fiberpoints[:, 1])
    fiberpoints[:, 0] += -np.amin(fiberpoints[:, 0])
    # threshdat.threshold*=256/threshdat.threshold.max()
    weights2d = threshdat.query(
        'fiber_diam==3 and sample==@samp2d').sort_values('master_fiber_index').threshold
    weights3d = threshdat.query(
        'fiber_diam==3 and sample==@samp3d').sort_values('master_fiber_index').threshold
    fiber2d = np.concatenate(
        [np.array(weights2d)[:, None], fiberpoints], axis=1)
    fiber3d = np.concatenate(
        [np.array(weights3d)[:, None], fiberpoints], axis=1)
    fiber2d = fiber2d.astype(np.float32)
    fiber3d = fiber3d.astype(np.float32)
    dist, _, flow = EMD(fiber2d, fiber3d, cv2.DIST_L2)
    print(f'dist for {samp3d}:', dist)

    def plot_flow(sig1, sig2, flow, arrow_width_scale=3):
        """Plots the flow computed by cv2.EMD.

        The source images are retrieved from the signatures and
        plotted in a combined image, with the first image in the
        red channel and the second in the green. Arrows are
        overplotted to show moved earth, with arrow thickness
        indicating the amount of moved earth.
        """

        img1 = sig_to_img(sig1)
        img2 = sig_to_img(sig2)
        combined = np.dstack((img1, img2, 0*img2))
        # RGB values should be between 0 and 1
        combined /= combined.max()
        print('Red channel is "before"; green channel is "after"; yellow means "unchanged"')
        plt.imshow(combined)

        flows = np.transpose(np.nonzero(flow))
        for src, dest in flows:
            # Skip the pixel value in the first element, grab the
            # coordinates. It'll be useful later to transpose x/y.
            start = sig1[src, 1:][::-1]
            end = sig2[dest, 1:][::-1]
            if np.all(start == end):
                # Unmoved earth shows up as a "flow" from a pixel
                # to that same exact pixel---don't plot mini arrows
                # for those pixels
                continue

            # Add a random shift to arrow positions to reduce overlap.
            shift = np.random.random(1) * .3 - .15
            start = start + shift
            end = end + shift

            mag = flow[src, dest] * arrow_width_scale
            plt.quiver(*start, *(end - start), angles='xy',
                       scale_units='xy', scale=1, color='white',
                       edgecolor='white', linewidth=mag/3,
                       width=mag, units='dots',
                       headlength=5,
                       headwidth=3,
                       headaxislength=4.5)

        plt.title("Earth moved from img1 to img2")

    def sig_to_img(sig):
        """Convert a signature back to a 2D image."""
        intsig = sig.astype(int)
        img = np.empty(
            (intsig[:, 1].max()+1, intsig[:, 2].max()+1), dtype=float)
        for i in range(sig.shape[0]):
            img[intsig[i, 1], intsig[i, 2]] = sig[i, 0]
        return img
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.imshow(sig_to_img(fiber2d), vmax=np.amax(
        fiber2d), vmin=np.amin(fiber2d), cmap='gray')
    ax1.set_title("Starting Distribution")
    ax2.imshow(sig_to_img(fiber3d))
    ax2.set_title("Ending Distribution")
    plt.show()

    plot_flow(fiber2d, fiber3d, flow)
    # %%
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Create a figure and axis
    plt.figure()

    # Create a Seaborn histplot with rug
    sns.histplot(fiber2d[:, 0], label="2DEM",
              fill=False,color='blue')
    sns.rugplot(fiber2d[:, 0],color='blue')
    sns.histplot(fiber3d[:, 0], label="3DM",  fill=False,color='orange')
    sns.rugplot(fiber3d[:, 0],color='orange')

    # Add a legend
    plt.legend()

    # Set labels for x and y axes
    plt.xlabel('Threshold (mA)')
    plt.ylabel('Density')

    # Show the plot
    plt.show()
