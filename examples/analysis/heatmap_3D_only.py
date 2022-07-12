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
# plt.rcParams.update({'font.size': 22})

# set default fig size
# plt.rcParams['figure.figsize'] = [16.8/3, 10.14*2 * 0.9]

samples = [2792]
samp3 = 2793
model = 0
sim = 33
sample_labels = ['Rostral\nContact','Caudal\n(Cathodic)\nContact']
nerve_label = '4R'
outpath = r'D:\ASCENT\ascent\out\analysis'
bigcomp = {0:'Anodic Leading',1:'Cathodic Leading'}

print("Note: assumed cathodic sample is the last one")
    
for i in range(len(samples)):
    # initialize and run Querys
    q = Query({
        'partial_matches': True,
        'include_downstream': True,
        'indices': {
            'sample': samples[i],
            'model': [model],
            'sim': [sim]
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
                title_toggle=True,
                title_override = 'Activation thresholds (mA) for mri 279',
                track_colormap_bounds=True,
            #    track_colormap_bounds_offset_ratio=0.0,
                # colomap_bounds_override=[[2.9,5.6],[2.9,5.6]],
            #    subplot_title_toggle=False,
                colorbar_text_size_override=30,
                # add_colorbar = False,
                # comp_sim = True,
                # override_axes =axs[i][1],
                dotsize = 15,
                # spec_nsim=2,
                show_orientation_point = True,
                thresh_source_sample = [samp3],
                #    tick_bounds=True,
                save_path = 'out/analysis'
           )
