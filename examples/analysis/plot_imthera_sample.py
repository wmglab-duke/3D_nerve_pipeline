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

import matplotlib.pyplot as plt
import seaborn as sns
from src.core import Sample
from src.core.query import Query
from src.utils import Object

criteria = {
    'partial_matches': True,
    'include_downstream': False,
    'indices': {'sample': [2515], 'model': None, 'sim': None},
}


q = Query(criteria)
q.run()

results = q.summary()

sample_index = results['samples'][0]['index']

fig, ax = plt.subplots(1, 1)
item: Sample = q.get_object(Object.SAMPLE, [results['samples'][0]['index']])
slide = item.slides[0]
slide.plot(
    fix_aspect_ratio=True,
    final=False,
    ax=ax,
    inner_index_labels=True,
    # scalebar=True,
    # scalebar_length=100,
    # scalebar_units='μm',
    fascicle_colors=sns.color_palette('Set1')[:4],
    inner_format='w-',
)
plt.xlabel('\u03bcm')
plt.ylabel('\u03bcm')

# add a circular band around the nerve with radius 1500 um from 0 to 270 degrees
theta = np.linspace(0, np.radians(325), 100)
r = 1700
x = r * np.cos(theta)
y = r * np.sin(theta)
plt.plot(x, y, 'k', linewidth=10, label='cuff')
# add 6 evenly spaced stars from 30 to 300 degrees, add numbers 0-5 to the stars
theta = np.linspace(np.radians(50), np.radians(325 - 50), 6)
r = 1550
x = r * np.cos(theta)
y = r * np.sin(theta)
plt.plot(x, y, 'ro', markersize=15, color='r', label='contacts')
for i in range(6):
    plt.text(x[i], y[i], str(i), fontsize=10, color='k', ha='center', va='center')
plt.legend(loc='lower right')


plt.show()
fname = str(sample_index)
fmt = 'svg'

dest = os.path.join('data', 'tmp', 'samples')
if not os.path.exists(dest):
    os.mkdir(dest)

fig.savefig(os.path.join(dest, f'{fname}.{fmt}'), format=fmt, dpi=1200)
