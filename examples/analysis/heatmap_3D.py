#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

import os
import sys

sys.path.append(os.path.sep.join([os.getcwd(), '']))

os.chdir('D:/ascent/ascent')

import matplotlib.pyplot as plt

from src.core.query import Query

# set default fig size
# plt.rcParams['figure.figsize'] = [16.8/3, 10.14*2 * 0.9]

samples = [250, 252]
samp3 = 253
model = 0
sim = 33
sample_labels = ['Rostral\nContact', 'Caudal\n(Cathodic)\nContact']
nerve_label = '2L'
outpath = r'D:\ASCENT\ascent\out\analysis'
bigcomp = {0: 'Anodic Leading', 1: 'Cathodic Leading'}

print("Note: assumed cathodic sample is the last one")

#%% Generates regular heatmaps with all nsims for each sample
allbound = {}

for i in range(len(samples)):
    # initialize and run Querys
    q = Query(
        {
            'partial_matches': True,
            'include_downstream': True,
            'indices': {'sample': samples[i], 'model': [model], 'sim': [sim]},
        }
    ).run()

    colormap_bounds_override = None

    useless, uselesstoo, bounds, labels = q.heatmaps(
        plot=False,
        # save_path='out/analysis',
        plot_mode='fibers',
        #    rows_override=6,
        colorbar_aspect=5,
        colormap_str='viridis',
        tick_count=4,
        reverse_colormap=True,
        title_toggle=True,
        title_override='2D extrusion thresholds: {} {}'.format(nerve_label, sample_labels[i].replace('\n', ' ')),
        track_colormap_bounds=True,
        #    track_colormap_bounds_offset_ratio=0.0,
        # colomap_bounds_override=[[2.9,5.6],[2.9,5.6]],
        #    subplot_title_toggle=False,
        colorbar_text_size_override=30,
        # add_colorbar = False,
        # override_axes = axs[i][0],
        # spec_nsim=2,
        dotsize=15,
        show_orientation_point=True,
        save_path='out/analysis'
        #    tick_bounds=True
    )

    useless, uselesstoo, bounds3d, labels = q.heatmaps(
        plot=False,
        # save_path='out/analysis',
        plot_mode='fibers',
        #    rows_override=6,
        colorbar_aspect=5,
        colormap_str='viridis',
        tick_count=4,
        reverse_colormap=True,
        title_toggle=True,
        title_override='3D model thresholds: {} {}'.format(nerve_label, sample_labels[i].replace('\n', ' ')),
        track_colormap_bounds=True,
        #    track_colormap_bounds_offset_ratio=0.0,
        # colomap_bounds_override=[[2.9,5.6],[2.9,5.6]],
        #    subplot_title_toggle=False,
        colorbar_text_size_override=30,
        # add_colorbar = False,
        # comp_sim = True,
        # override_axes =axs[i][1],
        dotsize=15,
        # spec_nsim=2,
        show_orientation_point=True,
        thresh_source_sample=[samp3],
        #    tick_bounds=True,
        save_path='out/analysis',
    )

    allbound[i] = []
    allbound[i].append(bounds)
    allbound[i].append(bounds3d)

#%% generates the paired heatmaps for all nsim and only cathodic contact
import numpy as np

for sample_index, sample_name in bigcomp.items():
    finalbound = np.array(allbound[sample_index])
    maxbound = np.amax(finalbound, axis=0)
    minbound = np.amin(finalbound, axis=0)
    override = np.zeros(maxbound.shape)
    override[:, 0] = minbound[:, 0]
    override[:, 1] = maxbound[:, 1]

    # need to: loop through and generate the 2x2 comparison for each nsim, then a whole new one
    # which does every nsim for the caudal contact for both
    fig, axs = plt.subplots(override.shape[0], 3, figsize=(10, 20), gridspec_kw={'width_ratios': [6, 6, 1]})

    # initialize and run Querys
    q = Query(
        {
            'partial_matches': True,
            'include_downstream': True,
            'indices': {'sample': samples[sample_index], 'model': [model], 'sim': [sim]},
        }
    ).run()

    q.heatmaps(
        plot=False,
        # save_path='out/analysis',
        plot_mode='fibers',
        #    rows_override=6,
        colorbar_aspect=5,
        colormap_str='viridis',
        tick_count=4,
        reverse_colormap=True,
        title_toggle=False,
        track_colormap_bounds=True,
        #    track_colormap_bounds_offset_ratio=0.0,
        colomap_bounds_override=override,
        subplot_title_toggle=False,
        colorbar_text_size_override=30,
        add_colorbar=False,
        override_axes=[x[0] for x in axs],
        # spec_nsim=2,
        dotsize=10,
        show_orientation_point=True
        #    tick_bounds=True
    )

    q.heatmaps(
        plot=False,
        # save_path='out/analysis',
        plot_mode='fibers',
        #    rows_override=6,
        colorbar_aspect=5,
        colormap_str='viridis',
        tick_count=4,
        reverse_colormap=True,
        title_toggle=False,
        track_colormap_bounds=True,
        #    track_colormap_bounds_offset_ratio=0.0,
        colomap_bounds_override=override,
        subplot_title_toggle=False,
        colorbar_text_size_override=30,
        add_colorbar=True,
        # comp_sim = True,
        cbar_axs=[x[2] for x in axs],
        override_axes=[x[1] for x in axs],
        dotsize=10,
        # spec_nsim=2,
        show_orientation_point=True,
        thresh_source_sample=[samp3]
        #    tick_bounds=True
    )

    for nsim, ax in enumerate(axs):
        this = ax[0]
        this.axis('on')
        this.set_ylabel(labels[nsim].split(':')[1], fontsize=20, rotation=0)
        this.axes.xaxis.set_visible(False)
        this.spines['top'].set_visible(False)
        this.spines['right'].set_visible(False)
        this.spines['bottom'].set_visible(False)
        this.spines['left'].set_visible(False)
        this.get_yaxis().set_ticks([])
    axs[0][2].set_title('mA', fontsize=30)
    axs[0][0].set_title('2D', fontsize=30)
    axs[0][1].set_title('3D', fontsize=30)
    plt.subplots_adjust(left=0.15)
    fig.text(0.05, 0.55, u"fiber diameter (\u03BCm)", ha="center", va="center", rotation=90, fontsize=35)
    fig.suptitle(f'Activation thresholds for {nerve_label}\n{sample_name} Contact', fontsize=40)
    fig.savefig(f'out/analysis/{samples[sample_index]}-{samp3}_{model}_{sim}', dpi=400, bbox_inches="tight")

#%% generates the individual comparison between all samples for each nsim

for i in range(override.shape[0]):
    fig, axs = plt.subplots(len(samples), 2, figsize=(25, 20))
    fig.suptitle(f'Activation thresholds for {nerve_label}, {labels[i]}', fontsize=60)
    for s in range(len(samples)):

        # initialize and run Querys
        q = Query(
            {
                'partial_matches': True,
                'include_downstream': True,
                'indices': {'sample': samples[s], 'model': [model], 'sim': [sim]},
            }
        ).run()

        q.heatmaps(
            plot=False,
            # save_path='out/analysis',
            plot_mode='fibers',
            #    rows_override=6,
            colorbar_aspect=5,
            colormap_str='viridis',
            tick_count=4,
            reverse_colormap=True,
            title_toggle=False,
            # track_colormap_bounds=True,
            #    track_colormap_bounds_offset_ratio=0.0,
            colomap_bounds_override=override,
            subplot_title_toggle=False,
            colorbar_text_size_override=30,
            add_colorbar=False,
            override_axes=[axs[s][0]],
            spec_nsim=i,
            dotsize=20,
            show_orientation_point=True
            #    tick_bounds=True
        )

        q.heatmaps(
            plot=False,
            # save_path='out/analysis',
            plot_mode='fibers',
            #    rows_override=6,
            colorbar_aspect=5,
            colormap_str='viridis',
            tick_count=4,
            reverse_colormap=True,
            title_toggle=False,
            # track_colormap_bounds=True,
            #    track_colormap_bounds_offset_ratio=0.0,
            colomap_bounds_override=override,
            subplot_title_toggle=False,
            colorbar_text_size_override=30,
            add_colorbar=False,
            # comp_sim = True,
            # cbar_axs = [x[2] for x in axs],
            override_axes=[axs[s][1]],
            dotsize=20,
            spec_nsim=i,
            show_orientation_point=True,
            thresh_source_sample=[samp3]
            #    tick_bounds=True
        )
        [ax.axis('on') for axes in axs for ax in axes]
        # TODO: make the second column 3d saample over and over again
        # TODO: lop over nsims and print params in title
        axs[0][0].set_title('2D extrusion model', fontsize=50)
        axs[0][1].set_title('Full 3D model', fontsize=50)
        for sa in range(len(samples)):
            axs[sa][0].set_ylabel(sample_labels[sa], fontsize=50, rotation=0, labelpad=100)
        for ax in axs:
            for a in ax:
                a.axes.xaxis.set_visible(False)
                a.spines['top'].set_visible(False)
                a.spines['right'].set_visible(False)
                a.spines['bottom'].set_visible(False)
                a.spines['left'].set_visible(False)
                a.get_yaxis().set_ticks([])
        plt.subplots_adjust(wspace=0, hspace=0)
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.8, 0.15, 0.05, 0.7])
        cmap = plt.get_cmap('viridis')
        cmap = cmap.reversed()
        mp = plt.cm.ScalarMappable(cmap=cmap)
        mp.set_clim(override[i, 0], override[i, 1])

        cb = fig.colorbar(
            mappable=mp,
            orientation='vertical',
            aspect=20,
            format='%0.2f',
            cax=cbar_ax,
        )
        cb.ax.tick_params(labelsize=50)
        cb.ax.tick_params(size=15)
        cb.ax.tick_params(width=5)

        cbar_ax.set_title(r'mA', fontsize=50)
        colorbar_text_size_override = 100
        # colorbar font size
        # if colorbar_text_size_override is not None:
        #     cb.ax.tick_params(labelsize=colorbar_text_size_override if (
        #             colorbar_text_size_override is not None) else 25)
        # cb.update_ticks()
        fig.subplots_adjust(left=0.15)
        fig.subplots_adjust(top=0.9)
        fig.savefig(f'out/analysis/{samples[-1]}-{samp3}_{model}_{sim}_{i}', dpi=400, bbox_inches="tight")
