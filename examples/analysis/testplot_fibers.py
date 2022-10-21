# -*- coding: utf-8 -*-
"""Created on Wed Sep  7 19:58:14 2022.

@author: Daniel
"""
# TODO: make this into a video animation
import os

import numpy as np
import pickle
from matplotlib import pyplot as plt

fibers = []
fiberpath = r'C:\Users\Daniel\Documents\ASCENT\threed\input\2LDS5\samples\3D\models\0\sims\3\3D_fiberset'
slidespath = r'C:\Users\Daniel\Documents\ASCENT\threed\input\slides\2Lslides.obj'
# load pickled slidelist
with open(slidespath, 'rb') as f:
    slidelist = pickle.load(f)

# load each fiber file and append to list
for file in os.listdir(fiberpath):
    print(file)
    if file.endswith('.dat'):
        fibers.append(np.loadtxt(os.path.join(fiberpath, file), skiprows=1))

#%%

# loop through each slide, and calculate z position
for i, slide in enumerate(slidelist):
    slide.scale(0.5)
    slide.scale(1.2)
    zpos = i * 20  # 20 microns per slice
    # get list of x,y coordinates for each fiber at this z position
    xs = []
    ys = []
    for fiber in fibers:
        # find the index of the closest z value
        idx = (np.abs(fiber[:, 2] - zpos)).argmin()
        # get the x,y coordinates at that index
        x = fiber[idx, 0]
        y = fiber[idx, 1]
        # append to slide list
        xs.append(x)
        ys.append(y)
    # plot the slide and all fiber points
    plt.figure()
    plt.plot(xs, ys, '.')
    plt.title(f'Slide {i}-zpos{zpos}')
    slide.plot()
    plt.show()
