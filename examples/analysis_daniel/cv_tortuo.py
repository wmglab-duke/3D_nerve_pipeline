#!/usr/bin/env python3
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

os.chdir('../..')

from src.core.query import Query

tortuosities = [1, 1.01, 1.05, 1.1]
total_distance = 2  # meters
actual_distance = [total_distance / t for t in tortuosities]
diams = {10: 55.6, 11.5: 63.7, 12.8: 71.2, 14: 78.4, 15: 85.5, 16: 91.6}
eff = {t: {key: value / t for key, value in diams.items()} for t in tortuosities}

df = pd.DataFrame({'tortuosity': tortuosities, 'distance': actual_distance})

# sns.pointplot(
#     data=df,
#     y='distance',
#     x='tortuosity'
#     )
effdat = (
    pd.DataFrame(eff)
    .melt(ignore_index=False)
    .reset_index()
    .rename(columns={'index': 'diam', 'variable': 'tortuosity', 'value': 'cv'})
)

sns.lineplot(data=effdat, x='diam', y='cv', hue='tortuosity')

dats = []
for sample in [253, 273, 373, 573, 653, 673]:
    q3 = Query(
        {
            'partial_matches': True,
            'include_downstream': True,
            'indices': {'sample': [sample], 'model': [0], 'sim': [3]},
        }
    ).run()
    dats.append(q3.data(source_sample=sample - 1, tortuosity=True))

dat3d = pd.concat(dats).reset_index(drop=True)
plt.figure()
sns.histplot(data=dat3d, x='tortuosity')

plt.figure()
sns.lmplot(
    data=dat3d.query('nsim in [0,5]'),
    facet_kws={'sharex': False, 'sharey': False},
    col='nsim',
    hue='sample',
    y='threshold',
    x='tortuosity',
)
