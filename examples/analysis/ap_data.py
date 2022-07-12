#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

import os
import sys

sys.path.append(os.path.sep.join([os.getcwd(), '']))
os.chdir('D:/ascent/ascent')

import numpy as np

import matplotlib.pyplot as plt
from src.core.query import Query
sys.path.append(os.getcwd()+'/subrepos/nd_line')
from nd_line.nd_line import nd_line
# os.chdir('D:/ascent/fresh')

# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14*2]) / 2)

samp3d = 453

samples =  [452]

models = [0]

sims = [33]

q = Query({
    'partial_matches': False,
    'include_downstream': True,
    'indices': {
        'sample': samples,
        'model': models,
        'sim': sims
    }
}).run()        
        
dat2d = q.ap_data(
    delta_V=60,
    absolute_voltage=False,
    delete_vmtime = True)

dat3d = q.ap_data(
    delta_V=60,
    absolute_voltage=False,
    delete_vmtime = True,
    sample_override = samp3d)

