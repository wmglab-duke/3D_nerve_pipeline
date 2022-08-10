#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

import os
import sys

sys.path.append(r'D:\ASCENT\ascent')
os.chdir(r'D:\ASCENT/ascent')

sys.path.append(os.path.sep.join([os.getcwd(), '']))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb

from src.core.query import Query

samples3d = [2, 3]

models = [0]

sims = [33]

dats = []
for sample in samples3d:
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [sample], 'model': models, 'sim': sims},
        }
    ).run()

    dats.append(q.threshdat3d(meanify=False))
data = pd.concat(dats)
sb.set_theme(style="whitegrid")
g = sb.catplot(
    x="sample",
    y='threshold',
    hue='sim',
    col="nsim",
    data=data,
    kind="strip",
    height=5,
    aspect=0.4,
    linewidth=1,
    sharey=False,
)
plt.gcf().savefig('out/analysis/alsocvss.png', dpi=400)
# pd.concat(dats).to_csv('out/analysis/threshes.csv',index=False)
