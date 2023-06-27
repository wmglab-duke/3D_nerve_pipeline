#!/usr/bin/env python3.7

"""Plot a sample.

The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent.

RUN THIS FROM REPOSITORY ROOT
"""

import os

import matplotlib.pyplot as plt

os.chdir('../..')

from src.core import Sample
from src.core.query import Query
from src.utils import Object

criteria = {
    'partial_matches': True,
    'include_downstream': False,
    'indices': {'sample': [672, 6729, 6721], 'model': None, 'sim': None},
}


q = Query(criteria)
q.run()

results = q.summary()
for sampledata in results['samples']:
    sample_index = sampledata['index']

    fig, ax = plt.subplots(1, 1)
    item: Sample = q.get_object(Object.SAMPLE, [sample_index])
    slide = item.slides[0]
    slide.plot(
        fix_aspect_ratio=True,
        final=False,
        ax=ax,
        inner_index_labels=False,
        # scalebar=True,
        # scalebar_length=100,
        # scalebar_units='μm',
    )
    ax.axis('off')
    plt.xlabel('\u03bcm')
    plt.ylabel('\u03bcm')
    plt.show()

    fname = str(sample_index)
    fmt = 'svg'
