"""Created on Wed Sep  7 19:58:14 2022.

@author: Daniel
"""
# TODO: make this into a video animation
import os
import pickle

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

threed = False
os.chdir('../..')
threshdata = pd.read_csv('thresh_unmatched_sim3_og.csv').query('nsim==5')
zpos = [29150, 29050, 28950]
for sample, samplenum, samp2d in zip(["6L"], [653], [652]):
    fibers = []
    fiberpath = os.path.join(os.getcwd(), fr'samples\{samplenum}\models\0\sims\3\3D_fiberset')
    slidespath = os.path.join(os.getcwd(), fr'input\slides\{sample}slides.obj')
    # load pickled slidelist
    with open(slidespath, 'rb') as f:
        slidelist = pickle.load(f)

    mfis = []
    # load each fiber file and append to list
    for file in os.listdir(fiberpath):
        if file.endswith('.dat'):
            fibers.append(np.loadtxt(os.path.join(fiberpath, file), skiprows=1))
            mfis.append(int(file.replace('.dat', '')))

    # get fiber thresholds
    thisdat = threshdata.query(f'sample in [@samplenum, @samp2d]')
    fig, axs = plt.subplots(len(zpos), 1, figsize=(6, 10))
    # for each zpos, plot the slide from the zpos and plot the fiber thresholds
    for z in zpos:
        slidenum = int(z / 20)
        slide = slidelist[slidenum]

        thresh = []
        x = []
        y = []
        # plot fiber thresholds
        for mfi, fiber in zip(mfis, fibers):
            # find the closest point to the zpos
            zloc = np.argmin(np.abs(fiber[:, 2] - z))
            # get the xy coords of the fiber at that zpos
            xy = fiber[zloc, :2]
            x.append(xy[0])
            y.append(-xy[1])
            # get the threshold for that fiber
            threshes = thisdat.query(f'master_fiber_index=={mfi} and sample=={samplenum if threed else samp2d}')[
                'threshold'
            ].values
            assert len(threshes) == 1
            thresh.append(threshes[0])

        # now do same as above but with vertical subplots
        axs[zpos.index(z)].scatter(x, y, c=thresh, cmap='viridis_r')
        slide.plot(final=False, inner_format='k-', ax=axs[zpos.index(z)])
        # turn axis off
        axs[zpos.index(z)].set_ylabel(f'z={z}')
        # turn off ticks and axis box
        axs[zpos.index(z)].set_xticks([])
        axs[zpos.index(z)].set_yticks([])
        axs[zpos.index(z)].spines['right'].set_visible(False)
        axs[zpos.index(z)].spines['left'].set_visible(False)
        axs[zpos.index(z)].spines['top'].set_visible(False)
        axs[zpos.index(z)].spines['bottom'].set_visible(False)
    fig.show()
    # add colorbar
    fig.subplots_adjust(right=0.85, hspace=0)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(axs[0].collections[0], cax=cbar_ax)
    # label colorbar
    cbar_ax.set_ylabel('Threshold (mA)')
    axs[0].set_title('3DM' if threed else '2DEM')
