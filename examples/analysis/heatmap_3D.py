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

samples = [670,672]
samp3 = 673
sample_labels = ['Rostral\nContact','Caudal\n(Cathodic)\nContact']
fig,axs = plt.subplots(len(samples),2,figsize=(25, 20))

for i in range(len(samples)):
    # initialize and run Querys
    q = Query({
        'partial_matches': True,
        'include_downstream': True,
        'indices': {
            'sample': samples[i],
            'model': [0],
            'sim': [3]
        }
    }).run()
    
    colormap_bounds_override = None

    q.heatmaps(plot=False,
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
                # colomap_bounds_override=[[2.9,5.6],[2.9,5.6]],
            #    subplot_title_toggle=False,
                colorbar_text_size_override=30,
                add_colorbar = False,
                override_axes = axs[i][0],
                spec_nsim=0,
                dotsize = 15,
                show_orientation_point = True
                #    tick_bounds=True
           )
    
    q.heatmaps(plot=False,
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
                # colomap_bounds_override=[[2.9,5.6],[2.9,5.6]],
            #    subplot_title_toggle=False,
                colorbar_text_size_override=30,
                add_colorbar = False,
                comp_sim = True,
                override_axes =axs[i][1],
                dotsize = 15,
                spec_nsim=0,
                show_orientation_point = True,
                thresh_source_sample = [samp3]
                #    tick_bounds=True
           )
    
# TODO: make the second column 3d saample over and over again
# TODO: lop over nsims and print params in title
axs[0][0].set_title('2D extrusion model',fontsize = 50)
axs[0][1].set_title('Full 3D model',fontsize = 50)
for i in range(len(samples)):
    axs[i][0].set_ylabel(sample_labels[i],fontsize = 50,rotation = 0,labelpad = 100)
for ax in axs:
    for a in ax:
        a.axes.xaxis.set_visible(False)
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.spines['bottom'].set_visible(False)
        a.spines['left'].set_visible(False)
        a.get_yaxis().set_ticks([])
plt.subplots_adjust(wspace = 0)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.8, 0.15, 0.05, 0.7])

mp = plt.cm.ScalarMappable()
mp.set_clim(2.9,5.6)

cb= fig.colorbar(
    mappable=mp,
    orientation='vertical',
    aspect=20,
    format='%0.2f',
    cax = cbar_ax,
)
cb.ax.tick_params(labelsize = 50)
cb.ax.tick_params(size = 15)
cb.ax.tick_params(width = 5)

cbar_ax.set_title(r'mA',fontsize = 50)
colorbar_text_size_override = 100
# colorbar font size
# if colorbar_text_size_override is not None:
#     cb.ax.tick_params(labelsize=colorbar_text_size_override if (
#             colorbar_text_size_override is not None) else 25)
# cb.update_ticks()
fig.subplots_adjust(left=.15)
fig.subplots_adjust(top=.95)
fig.savefig('out/analysis/3dhm.png',dpi=500)
