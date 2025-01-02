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
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from src.core import Sample
from src.core.query import Query
from src.utils import Object

lim = [-2000, 2000]
max_thk = 1000

# sns.set(font_scale=1,style='white')
for samplenum, r_cuff_in_pre_MCT in zip([25212, 37212, 57212, 67212], [1000, 1500, 1000, 1000]):
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
    # theta = np.linspace(0, np.radians(325), 100)
    # r = 100+slide.nerve.ecd()/2
    R_in_MCT = slide.nerve.ecd() / 2 + 150
    thetamax = np.radians((1 * (r_cuff_in_pre_MCT / R_in_MCT)) * 360)
    theta = np.linspace(0, thetamax, 100)

    x = R_in_MCT * np.cos(theta)
    y = R_in_MCT * np.sin(theta)
    plt.plot(x, y, 'k', linewidth=10, label='cuff')
    # add 6 evenly spaced stars from 30 to 300 degrees, add numbers 0-5 to the stars
    theta = np.linspace(np.radians(50), np.radians(325 - 50), 4)
    # add a circular band around the nerve with a certain arc length (e.g., 10 mm)

    for i in range(4):
        ang_contactcenter_pre_MCT = 90
        ang_cuffseam_pre_MCT = 45
        ang_contactcenter_MCT = ang_contactcenter_pre_MCT * (r_cuff_in_pre_MCT / R_in_MCT)
        ang_cuffseam_MCT = ang_cuffseam_pre_MCT * (r_cuff_in_pre_MCT / R_in_MCT)
        contactpos = np.radians(ang_cuffseam_MCT + i * ang_contactcenter_MCT)
        contactspan = np.radians(10)
        contactheta = np.linspace(contactpos - contactspan, contactpos + contactspan, 100)
        r = 100 + slide.nerve.ecd() / 2
        x = r * np.cos(contactheta)
        y = r * np.sin(contactheta)
        plt.plot(x, y, color='r', label='contacts' if i == 0 else '_', linewidth=5)
        txt = plt.text(x[int(len(x) / 2)], y[int(len(x) / 2)], str(i), fontsize=12, color='k', ha='center', va='center')
        txt.set_path_effects([PathEffects.withStroke(linewidth=1.5, foreground='w')])

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
