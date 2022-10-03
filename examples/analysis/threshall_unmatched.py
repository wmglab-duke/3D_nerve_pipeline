#!/usr/bin/env python3
import json
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.use('agg')
sys.path.append('.')

from src.core.query import Query

datas = []
model = 0
source_sim = 3
with open('examples/analysis/plotconfig.json') as f:
    config = json.load(f)
for simdex in config['sim_data'].keys():
    simint = int(simdex)
    for sample_data in config['sample_data']:
        samp3d = sample_data['index3d']
        print(sample_data)
        nerve_label = sample_data['name']
        samples2d = [x['index'] for x in sample_data['exsamples']]
        q = Query(
            {
                'partial_matches': False,
                'include_downstream': True,
                'indices': {'sample': samples2d, 'model': [model], 'sim': [simint]},
            }
        ).run()
        dat2d = q.data(tortuosity=True, peri_site=True, zpos=True, cuffspan=[27000, 29000], label='2L')
        dat2d['type'] = '2D'
        dat2d['contact'] = ''
        anodic = dat2d['sample'].astype(str).str.endswith('0')
        dat2d.loc[anodic, 'contact'] = 'anodic'
        cathodic = dat2d['sample'].astype(str).str.endswith('2')
        dat2d.loc[cathodic, 'contact'] = 'cathodic'
        datas.append(dat2d)
        q3 = Query(
            {
                'partial_matches': False,
                'include_downstream': True,
                'indices': {'sample': [samp3d], 'model': [model], 'sim': [simint]},
            }
        ).run()
        dat3d = q3.data(
            source_sample=samples2d[0],
            tortuosity=True,
            peri_site=True,
            zpos=True,
            cuffspan=[28000, 30000],
            source_sim=source_sim,
            label='2L',
        )
        dat3d['type'] = '3D'
        datas.append(dat3d)
data = pd.concat(datas)
data.to_csv('thresh_unmatched.csv', index=False)
