"""Created on Wed Sep  7 19:58:14 2022.

@author: Daniel
"""
# TODO: make this into a video animation
import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

os.chdir('../..')

# for sample, samplenum in zip(["2L", "2R", "3R", "5R", "6L", "6R"], [253, 273, 373, 573, 653, 673]):
for sample, samplenum in zip(["5Rdef"], [5735]):
    fibers = {}
    fiberpath = os.path.join(os.getcwd(), fr'samples\{samplenum}\models\0\sims\3\3D_fiberset')
    slidespath = os.path.join(os.getcwd(), fr'input\slides\{sample}slides.obj')

    # load pickled slidelist
    with open(slidespath, 'rb') as f:
        slidelist = pickle.load(f)

    # load each fiber file and append to list
    for file in os.listdir(fiberpath):
        print(file)
        # if file == '1.dat':# use this to watch a specific fiber
        if file.endswith('.dat'):
            fibers[file.replace('.dat', '')] = np.loadtxt(os.path.join(fiberpath, file), skiprows=1)

    # %%

    # loop through each slide, and calculate z position
    for i, slide in enumerate(slidelist):
        plt.figure()
        # slide.scale(0.5)
        zpos = i * 20  # 20 microns per slice
        if zpos < 25000 or zpos > 26000:
            continue
        if i % 5 != 0:
            continue
        # if zpos<27000 or zpos>29000: continue #use this to check only a specific range
        # get list of x,y coordinates for each fiber at this z position
        xs = []
        ys = []
        for mfi, fiber in fibers.items():
            # find the index of the closest z value
            idx = (np.abs(fiber[:, 2] - zpos)).argmin()
            # get the x,y coordinates at that index
            x = fiber[idx, 0]
            y = fiber[idx, 1]
            plt.scatter(x, -y, s=3, color='red')
            plt.text(x, -y, mfi, fontsize=5)
        # plot the slide and all fiber points
        # plt.scatter(xs, ys, s=3, color='red')
        plt.title(f'Slide {i}-zpos{zpos}')
        slide.plot()
        plt.show()
