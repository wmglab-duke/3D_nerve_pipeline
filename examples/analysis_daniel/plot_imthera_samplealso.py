#!/usr/bin/env python3.7

"""Plot a sample.

The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent.

RUN THIS FROM REPOSITORY ROOT
"""

import os

import numpy as np

os.chdir('../..')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import to_hex
from src.core import Sample
from src.core.query import Query
from src.utils import Object

lim = [-2000, 2000]
max_thk = 1000

# sns.set(font_scale=1,style='white')
for samplenum in [25212, 37212, 57212, 67212]:
    criteria = {
        'partial_matches': True,
        'include_downstream': False,
        'indices': {'sample': [samplenum], 'model': None, 'sim': None},
    }

    q = Query(criteria)
    q.run()

    results = q.summary()

    sample_index = results['samples'][0]['index']

    fig, ax = plt.subplots(1, 1)
    item: Sample = q.get_object(Object.SAMPLE, [results['samples'][0]['index']])
    slide = item.slides[0]
    cmap = cm.get_cmap('rainbow', len(slide.inners()))
    slide.plot(
        fix_aspect_ratio=True,
        final=False,
        ax=ax,
        inner_index_labels=False,
        # scalebar=True,
        # scalebar_length=100,
        # scalebar_units='μm',
        fascicle_colors=[to_hex(cmap(i)) for i in range(len(slide.inners()))],
        inner_format='w-',
        show_axis=False,
    )
    plt.xlabel('\u03bcm')
    plt.ylabel('\u03bcm')
    plt.ylim(lim)
    plt.xlim(lim)
    slide.add_scalebar(plt.gca(), 1, 'mm')
    # add a circular band around the nerve with radius 1500 um from 0 to 270 degrees
    theta = np.linspace(0, np.radians(325), 100)
    r = 100 + slide.nerve.ecd() / 2
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    plt.plot(x, y, 'k', linewidth=5, label='cuff')
    # add 6 evenly spaced stars from 30 to 300 degrees, add numbers 0-5 to the stars
    theta = np.linspace(np.radians(50), np.radians(325 - 50), 4)
    r = 100 + slide.nerve.ecd() / 2
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    plt.plot(x, y, 'ro', markersize=15, color='r', label='contacts')
    for i in range(4):
        plt.text(x[i], y[i], str(i), fontsize=10, color='k', ha='center', va='center')
    plt.legend(loc='upper right')

    # plot fascicle ECD
    ECD = {i: f.ecd() for i, f in enumerate(slide.inners())}
    colors = [to_hex(cmap(i)) for i in range(len(slide.inners()))]
    # plot sxcatterplot of ECD with x=0 for all
    fig, ax = plt.subplots(1, 1)
    for i, ecd in ECD.items():
        ax.scatter(0, ecd, color=colors[i], marker='s', s=100)
        ax.scatter(0, ecd, color='black', marker='_', s=100)
    plt.gca().set_aspect(0.0015)
    plt.xticks([], [])
    plt.ylabel('fascicle diameter (μm)')
    plt.ylim(0, max_thk)
    plt.show()
