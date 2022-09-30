#!/usr/bin/env python3
import json
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from src.core.query import Query

sys.path.insert(0, os.path.abspath('../../'))
os.chdir('../../')
from src.core.plotter import get_datamatch

samples2d = [250, 252]
samp3d = 253
model = 0
source_sim = 3
simdex = 3
nerve_label = '2L'
q = Query(
    {
        'partial_matches': False,
        'include_downstream': True,
        'indices': {'sample': samples2d, 'model': [model], 'sim': [simdex]},
    }
).run()
dat2d = q.data(tortuosity=True, peri_site=True, zpos=True, cuffspan=[27000, 29000], label='2L')
q3 = Query(
    {
        'partial_matches': False,
        'include_downstream': True,
        'indices': {'sample': [samp3d], 'model': [model], 'sim': [simdex]},
    }
).run()
dat3d = q3.data(
    source_sample=samples2d[0],
    tortuosity=True,
    peri_site=True,
    zpos=True,
    cuffspan=[27000, 29000],
    source_sim=source_sim,
    label='2L',
)
# TODO add note cuff cathodic anodic
