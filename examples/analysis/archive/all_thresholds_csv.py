#!/usr/bin/env python3
import json
import sys

import matplotlib
import pandas as pd

matplotlib.use('agg')
sys.path.append('.')

from src.core.plotter import get_datamatch

datas = []
model = 0
with open('examples/analysis/plotconfig.json') as f:
    config = json.load(f)
for simdex in config['sim_data'].keys():
    simint = int(simdex)
    for sample_data in config['sample_data']:
        print(simdex, sample_data)
        samp3d = sample_data['index3d']
        nerve_label = sample_data['name']
        samples2d = [x['index'] for x in sample_data['exsamples']]
        threshdat = get_datamatch(
            samples2d, samp3d, model, simint, nerve_label, source_sim=3, tortuosity=True, cuffspan=[27000, 29000]
        )
        datas.append(threshdat)
data = pd.concat(datas)
data.to_csv('thresh.csv', index=False)
