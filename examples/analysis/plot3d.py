#!/usr/bin/env python3.7

"""The copyrights of this software are owned by Duke University.

Please refer to the LICENSE.txt and README.txt files for licensing
instructions. The source code can be found on the following GitHub
repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

import json
import os
import sys

sys.path.append('.')

from src.core.plotter import (
    ap_plot,
    plot_colorjoint,
    plot_colorthresh,
    plot_colorthreshstrip,
    plot_correlation,
    plot_dose_response,
    plot_oneone,
)

model = 0
cuff_contacts = [21, 29]

with open('examples/analysis/plotconfig.json') as f:
    config = json.load(f)

for simdex in config['sim_data'].keys():
    simint = int(simdex)
    os.makedirs(f'out/analysis/{simdex}', exist_ok=True)
    for sample_data in config['sample_data']:
        samp3d = sample_data['index3d']
        nerve_label = sample_data['name']
        samples2d = [x['index'] for x in sample_data['exsamples']]
        plot_correlation(samples2d, samp3d, model, simint, nerve_label)
        plot_oneone(samples2d, samp3d, model, simint, nerve_label)
        plot_dose_response(samples2d, samp3d, model, simint, nerve_label)
        for samp2d in samples2d:
            plot_colorthreshstrip(samp2d, samp3d, model, simint, nerve_label)
            plot_colorthresh(samp2d, samp3d, model, simint, nerve_label)
            plot_colorjoint(samp2d, samp3d, model, simint, nerve_label)
            ap_plot(samp2d, samp3d, model, simint, cuff_contacts, source_sim=config['source_sim'])
