"""Created on Fri Jun 16 14:11:57 2023.

@author: dpm42
"""

import json
import os
import pickle

import numpy as np
from src.core import Sample
from src.utils import Config, SetupMode, WriteMode

for sample in os.listdir('samples'):
    if len(sample) < 3 or str(sample)[2] == '3' or sample.startswith('.'):
        continue
    # load sample object
    sample_obj = pickle.load(open(os.path.join('samples', sample, 'sample.obj'), 'rb'))
    refslide = sample_obj.slides[0]
    # create slide from input
    config = json.load(open(os.path.join('samples', sample, 'sample.json')))
    print('SHRINKAGE:', config['scale']['shrinkage'])
    run_config = json.load(open(os.path.join('config', 'user', 'runs', sample + '.json')))
    old_area = refslide.nerve.area()
    sample_obj.slides = []
    sample_obj.populate()
    new_area = sample_obj.slides[0].nerve.area()
    print(f'{sample}: {old_area} -> {new_area}')
    assert np.isclose(old_area, new_area, rtol=1e-3, atol=1e-3)
