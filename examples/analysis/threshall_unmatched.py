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
for simdex in ['3']:
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
        dat2d = q.data(tortuosity=True, peri_site=True, zpos=True, cuffspan=[28000, 30000], label=nerve_label)
        dat2d['type'] = '2D'
        dat2d['contact'] = ''
        anodic = dat2d['sample'].astype(str).str.endswith('0')
        dat2d.loc[anodic, 'contact'] = 'anodic'
        cathodic = dat2d['sample'].astype(str).str.endswith('2')
        dat2d.loc[cathodic, 'contact'] = 'cathodic'
        q3 = Query(
            {
                'partial_matches': False,
                'include_downstream': True,
                'indices': {'sample': [samp3d], 'model': [model], 'sim': [simint]},
            }
        ).run()
        dat3d = q3.data(
            source_sample=samples2d[0],
            # tortuosity=True,
            # peri_site=True,
            # zpos=True,
            # cuffspan=[28000, 30000],
            source_sim=source_sim,
            label=nerve_label,
        )
        dat3d['type'] = '3D'
        # for each nsim within the sim, use the fiber_diam and pulse_width from the nsim key
        # each nsim has a different pulse width and fiber diameter
        nsim_key = config['sim_data'][simdex]['nsim_key']
        for nsim in nsim_key:
            pulse_width = nsim_key[nsim]['pulse_width']
            fiber_diam = nsim_key[nsim]['fiber_diam']
            # find where dat2d and dat3d sim and nsim match the current sim and nsim
            # then add the pulse width and fiber diameter to the data
            dat2d.loc[(dat2d['sim'] == simint) & (dat2d['nsim'] == int(nsim)), 'pulse_width'] = pulse_width
            dat2d.loc[(dat2d['sim'] == simint) & (dat2d['nsim'] == int(nsim)), 'fiber_diam'] = fiber_diam
            dat3d.loc[(dat3d['sim'] == simint) & (dat3d['nsim'] == int(nsim)), 'pulse_width'] = pulse_width
            dat3d.loc[(dat3d['sim'] == simint) & (dat3d['nsim'] == int(nsim)), 'fiber_diam'] = fiber_diam
        datas.append(dat2d)
        datas.append(dat3d)
data = pd.concat(datas)
data.to_csv('thresh_unmatched.csv', index=False)
