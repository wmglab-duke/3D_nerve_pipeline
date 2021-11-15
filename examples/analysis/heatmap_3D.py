#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

# ASSUMPTIONS: (old, moved to Query.heatmaps)
#   1) 1:1 inner:outer for all fascicles
#   2) Single slide for each sample (0_0)
#   3) Single fiber per inner
# TODO: Change above assumptions in later iteration? (highest priority is probably assumption 3)

import os
import sys

sys.path.append(os.path.sep.join([os.getcwd(), '']))

os.chdir('D:/ascent/fresh')

import matplotlib.pyplot as plt
from src.core.query import Query

# set default fig size
plt.rcParams['figure.figsize'] = [16.8/3, 10.14*2 * 0.9]
fig,axs = plt.subplots(3,2,figsize=(25, 25))

samples = [3070,3150,3230]
for i in range(3):
    # initialize and run Querys
    q = Query({
        'partial_matches': True,
        'include_downstream': True,
        'indices': {
            'sample': samples[i],
            'model': [0],
            'sim': [3000]
        }
    }).run()
    
    # NOTE: these values were copied from the output of heatmaps(), setting the track_colormap_bounds flag True
    colormap_bounds_override = None
# colormap_bounds_override = [
#     [
#         (0.322187, 1.51728),
#         (0.110607, 0.64657),
#         (0.03439, 0.133777),
#         (0.026082, 0.091096),
#         (0.023491, 0.077376),
#         (0.020899, 0.06945),
#     ],
#     [
#         (0.568523, 1.4246),
#         (0.177679, 0.624619),
#         (0.049023, 0.161826),
#         (0.036371, 0.115485),
#         (0.031646, 0.099632),
#         (0.027606, 0.089571),
#     ],
#     [
#         (0.658765, 1.619717),
#         (0.206337, 0.727056),
#         (0.059389, 0.19963),
#         (0.044908, 0.15024),
#         (0.039725, 0.131948),
#         (0.034542, 0.118534),
#     ],
#     [
#         (0.946562, 1.4246),
#         (0.280725, 0.62218),
#         (0.068535, 0.152679),
#         (0.048566, 0.107559),
#         (0.040335, 0.092315),
#         (0.035304, 0.082864),
#     ]
# ]

# colormap_bounds_override = [
#     [
#         (0.35440570000000005, 1.365552),
#         (0.1216677, 0.581913),
#         (0.037829, 0.12039930000000001),
#         (0.028690200000000003, 0.0819864),
#         (0.025840100000000005, 0.0696384),
#         (0.022988900000000003, 0.062505),
#     ],
#     [
#         (0.6253753000000001, 1.28214),
#         (0.1954469, 0.5621571000000001),
#         (0.0539253, 0.1456434),
#         (0.040008100000000005, 0.1039365),
#         (0.034810600000000004, 0.0896688),
#         (0.0303666, 0.0806139),
#     ],
#     [
#         (0.7246415000000002, 1.4577453),
#         (0.2269707, 0.6543504),
#         (0.06532790000000001, 0.179667),
#         (0.0493988, 0.135216),
#         (0.04369750000000001, 0.11875320000000002),
#         (0.03799620000000001, 0.1066806),
#     ],
#     [
#         (1.0412182, 1.28214),
#         (0.3087975, 0.559962),
#         (0.07538850000000001, 0.1374111),
#         (0.0534226, 0.0968031),
#         (0.044368500000000005, 0.08308349999999999),
#         (0.038834400000000005, 0.0745776),
#     ]
# ]

# builds heatmaps


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
                colomap_bounds_override=[[2.9,5.6],[2.9,5.6]],
            #    subplot_title_toggle=False,
                colorbar_text_size_override=30,
                add_colorbar = False,
                comp_sim = True,
                override_axes =axs[i],
                dotsize = 15,
                show_orientation_point = True
            #    tick_bounds=True
           )

#
#                 # TODO: Finish building heatmap of polyfasc nerve (1 fiber/fasc)
#                 # also, look into adding documentation to Simulation (might be useful for above task too)

#plt.close('all')

axs[0][0].set_title('2D extrusion model',fontsize = 50)
axs[0][1].set_title('Full 3D model',fontsize = 50)
axs[0][0].set_ylabel('Rostral\nContact',fontsize = 50,rotation = 0,labelpad = 100)
axs[1][0].set_ylabel('Center',fontsize = 50,rotation = 0,labelpad = 100)
axs[2][0].set_ylabel('Caudal\nContact',fontsize = 50,rotation = 0,labelpad = 100)
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
