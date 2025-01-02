#!/usr/bin/env python3
import json
import os
import sys

import pandas as pd

sys.path.append('.')
import os

# os.chdir('../..')
from src.core.query import Query

model = 0
source_sim = 3
with open('examples/analysis/plotconfig_og_mono.json') as f:
    config = json.load(f)
for simdex in config['sim_data'].keys():
    simint = int(simdex)
    datas = []
    for sample_data in config['sample_data']:
        samp3d = sample_data['index3d']
        print(sample_data)
        nerve_label = sample_data['name']
        samples2d = sample_data['exsamples']
        q = Query(
            {
                'partial_matches': False,
                'include_downstream': True,
                'indices': {'sample': samples2d, 'model': [model], 'sim': [simint]},
            }
        ).run()
        dat2d = q.data(
            tortuosity=True,
            peri_site=True,
            zpos=True,
            cuffspan=[20000, 30000],
            label=nerve_label,
            efib_distance=True,
            oneten=True,
        )
        dat2d['type'] = '2D'
        dat2d['contact'] = ''
        anodic = dat2d['sample'].astype(str).str[2] == '0'
        dat2d.loc[anodic, 'contact'] = 'anodic'
        cathodic = dat2d['sample'].astype(str).str[2] == '2'
        dat2d.loc[cathodic, 'contact'] = 'cathodic'
        q3 = Query(
            {
                'partial_matches': False,
                'include_downstream': True,
                'indices': {'sample': [samp3d], 'model': [model], 'sim': [simint]},
            }
        ).run()
        source_samples = [x for x in samples2d if str(x)[2] == '2']  # use cathodic sample as source
        assert len(source_samples) == 1
        dat3d = q3.data(
            source_sample=source_samples[0],
            tortuosity=True,
            peri_site=True,
            zpos=True,
            cuffspan=[20000, 30000],
            source_sim=source_sim,
            label=nerve_label,
            efib_distance=True,
            oneten=True,
        )
        dat3d['type'] = '3D'
        # for each nsim within the sim, use the fiber_diam and pulse_width from the nsim key
        # each nsim has a different pulse width and fiber diameter
        nsim_key = config['sim_data'][simdex]['nsim_key']
        for nsim in nsim_key.keys():
            pulse_width = nsim_key[nsim]['pulse_width']
            fiber_diam = nsim_key[nsim]['fiber_diam']
            # find where dat2d and dat3d sim and nsim match the current sim and nsim
            # then add the pulse width and fiber diameter to the data
            dat2d.loc[(dat2d['sim'] == simint) & (dat2d['nsim'] == nsim), 'pulse_width'] = pulse_width
            dat2d.loc[(dat2d['sim'] == simint) & (dat2d['nsim'] == nsim), 'fiber_diam'] = fiber_diam
            dat3d.loc[(dat3d['sim'] == simint) & (dat3d['nsim'] == nsim), 'pulse_width'] = pulse_width
            dat3d.loc[(dat3d['sim'] == simint) & (dat3d['nsim'] == nsim), 'fiber_diam'] = fiber_diam
        datas.append(dat2d)
        datas.append(dat3d)
    data = pd.concat(datas)
    data.to_csv(f'thresh_unmatched_sim{simint}_og_mono.csv', index=False)
