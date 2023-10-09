"""Created on Wed Jan 19 12:55:26 2022.

@author: dpm42
"""
import os
import shutil

import numpy as np

os.chdir('D:/ascent/ascent/samples')

samples = [473, 453]

model = 0

source_sim = 3

dest_sim = 33

for sample in samples:
    sim_dir = f'{sample}/models/{model}/sims/{source_sim}'
    dest_dir = f'{sample}/models/{model}/sims/{dest_sim}'
    for name in ['ss_bases', 'ss_coords', 'ss_lengths']:
        shutil.copytree(os.path.join(sim_dir, name), os.path.join(dest_dir, name))
    for name in ['ss_bases/0', 'ss_bases/1']:
        for file in os.listdir(os.path.join(dest_dir, name)):
            data = np.loadtxt(os.path.join(dest_dir, name, file), skiprows=1)
            data = np.flip(data, axis=0)
            np.savetxt(os.path.join(dest_dir, name, file), data, header=str(len(data)), comments='')
