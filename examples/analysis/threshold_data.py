#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

import os
import sys
sys.path.append(r'D:\ASCENT\fresh')
os.chdir(r'D:\ASCENT/fresh')

sys.path.append(os.path.sep.join([os.getcwd(), '']))

import numpy as np

import matplotlib.pyplot as plt
from src.core.query import Query
import pandas as pd

threed = 19

samples = [19]

models = [451]

sims = [1819]

dats = []

q = Query({
    'partial_matches': False,
    'include_downstream': True,
    'indices': {
        'sample': samples,
        'model': models,
        'sim': sims
    }
}).run()

# builds heatmaps
# q.barcharts_compare_models(logscale=False,
#                            model_labels=['Model 0: Veltink Epineurium, \n              Veltink Perineurium',
#                                          'Model 1: Veltink Epineurium, \n              Goodall Perineurium',
#                                          'Model 2: Goodall Epineurium, \n              Veltink Perineurium',
#                                          'Model 3: Goodall Epineurium, \n              Goodall Perineurium']
#                            )
dat2d = q.threshdat(sl=False,meanify=False)

q = Query({
    'partial_matches': False,
    'include_downstream': True,
    'indices': {
        'sample': [threed],
        'model': models,
        'sim': sims
    }
}).run()

dat3d = q.threshdat3d(meanify = False)

dat2d['threshold'] = np.nan
for i in range(len(dat2d)):
    row = dat2d.iloc[i,:]
    thresh = dat3d[(dat3d["sample"]==row['sample']) & 
                       (dat3d["model"]==row['model']) &
                       (dat3d["sim"]==row['sim']) &
                       (dat3d["nsim"]==row['nsim']) &
                       (dat3d["index"]==row['index'])]['threshold']
    thresh=list(thresh)
    if len(thresh)!=1:sys.exit('issue here')
    dat2d.iloc[i,-1]=thresh[0]
if np.any(dat2d.threshold==np.nan):
    sys.exit('issue here too')


# pd.concat(dats).to_csv('out/analysis/threshes.csv',index=False)
