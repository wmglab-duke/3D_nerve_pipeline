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

criteria = {
    'partial_matches': True,
    'include_downstream': False,
    'indices': {'sample': [6729], 'model': None, 'sim': None},
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
    inner_index_labels=True,
    # scalebar=True,
    # scalebar_length=100,
    # scalebar_units='μm',
    fascicle_colors=[to_hex(cmap(i)) for i in range(len(slide.inners()))],
    inner_format='w-',
)
plt.xlabel('\u03bcm')
plt.ylabel('\u03bcm')

# load contact coordinates
dir = r'D:\threed_final\input\contact_coords\6RdefDS5'
for i in [1]:  # NOTE set this to only load caudal coords
    contact_coords = np.loadtxt(dir + r'\pcs' + str(i + 1) + '.txt', skiprows=8)[:, :2]
    plt.scatter(contact_coords[:, 0], contact_coords[:, 1], color='k', label='contacts' if i == 0 else '_')
    plt.legend()
plt.show()
fname = str(sample_index)
fmt = 'svg'

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
plt.ylim(0, 1000)
plt.show()
