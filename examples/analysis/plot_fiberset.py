#!/usr/bin/env python3.7

"""Generate a plot of fiber coordinates overlaid with a plot of the sample.

The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent.

RUN THIS FROM REPOSITORY ROOT
"""

import os

os.chdir('../..')

import matplotlib.pyplot as plt

from src.core import Sample, Simulation
from src.core.query import Query
from src.utils import Object

criteria = {
    'partial_matches': True,
    'include_downstream': True,
    'indices': {'sample': [572], 'model': [0], 'sim': [3]},
}


q = Query(criteria)
q.run()

results = q.summary()

sample_index = results['samples'][0]['index']
model_index = results['samples'][0]['models'][0]['index']
sim_index = results['samples'][0]['models'][0]['sims'][0]

sample: Sample = q.get_object(Object.SAMPLE, [results['samples'][0]['index']])
sim: Simulation = q.get_object(Object.SIMULATION, [sample_index, model_index, sim_index])

for fiberset_ind, fiberset in enumerate(sim.fibersets):
    if fiberset_ind != 0:
        continue
    slide = sample.slides[0]
    fig, ax = plt.subplots(1, 1)
    slide.plot(fix_aspect_ratio=True, final=False, ax=ax)
    fiberset.plot(
        ax=ax, indices=True, inner_specific_indices=False, scatter_kws={'s': 1}, annotate_kws=dict(ha='center', size=4)
    )

    plt.xlabel('\u03bcm')
    plt.ylabel('\u03bcm')

    fname = f'{str(sample_index)}_{str(model_index)}_{str(sim_index)}_{str(fiberset_ind)}'
    fmt = 'png'
    plt.gcf().savefig('testout.png', dpi=400)
    print(sim.indices_fib_to_n(0, 26))
