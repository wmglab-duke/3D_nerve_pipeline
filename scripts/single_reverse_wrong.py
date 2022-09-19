# -*- coding: utf-8 -*-
"""Created on Wed Jan 19 12:55:26 2022.

@author: dpm42
"""
import os

import numpy as np

path = r'D:\ASCENT\ascent\samples\473\models\0\sims\33\3D_fiberset'
for file in os.listdir(path):
    data = np.loadtxt(os.path.join(path, file), skiprows=1)
    data = np.flip(data, axis=0)
    np.savetxt(os.path.join(path, file), data, header=str(len(data)), comments='')
