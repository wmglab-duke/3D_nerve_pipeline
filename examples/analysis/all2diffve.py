import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

os.chdir('../../')
model = 0
sim = 3
nsim = 0

for sample, samp3d, nerve_label in zip([372, 3721], [373, 3731], ['3R', '3Rdef']):
    threshdat = pd.read_csv('thresh_unmatched_sim3_og.csv')
    threshdat = threshdat.query(f'sample=={sample} & nsim==0 & model=={model} & sim=={sim}')
    threshdat.reset_index(drop=True, inplace=True)
    sns.set_style('whitegrid')

    base_n_sim = os.path.join('samples', str(sample), 'models', str(model), 'sims', str(sim), 'n_sims')
    base_n_sim3d = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'n_sims')

    # make figure with 2x2 subplots, share y axis
    fig, axs = plt.subplots(2, 1, sharey=True, sharex=True)

    fig2, ax2 = plt.subplots()

    # loop through each fiber
    for i, row in threshdat.iterrows():
        inner = int(row['inner'])
        fiber = int(row['fiber'])
        inner3d = 0
        fiber3d = int(row['master_fiber_index'])
        pve2 = os.path.join(base_n_sim, str(nsim), 'data', 'inputs', f'inner{inner}_fiber{fiber}.dat')
        pve3 = os.path.join(base_n_sim3d, str(nsim), 'data', 'inputs', f'inner{inner3d}_fiber{fiber3d}.dat')
        v2 = np.loadtxt(pve2)
        v3 = np.loadtxt(pve3)
        v2 = v2[1::11]
        v3 = v3[1::11]
        der2 = np.diff(np.diff(v2))
        der3 = np.diff(np.diff(v3))
        der2 = der2 / np.max(np.abs(der2))
        der3 = der3 / np.max(np.abs(der3))
        der3 = der3[::-1]

        # plot second differential in the first subplot
        axs[0].plot(der2, alpha=0.1, color='r')
        # plot base Ve values in the second subplot
        # ax2.plot(v2/np.amax(np.abs(v2)), alpha=0.1, color='b')
        ax2.plot(v2, alpha=0.1, color='b')

        # plot second differential in the third subplot
        axs[1].plot(der3, alpha=0.1, color='r')
        # plot base Ve values in the fourth subplot
        # ax2.plot(v3[::-1]/np.amax(np.abs(v3)), alpha=0.1, color='k')
        ax2.plot(v3[::-1], alpha=0.1, color='k')

        print(i)

    # axs[0, 0].set_ylabel('Normalized 2nd Differential')
    # axs[1, 0].set_ylabel('Normalized 2nd Differential')
    # axs[0, 1].set_ylabel('Ve')
    # axs[1, 1].set_ylabel('Ve')
    # axs[1, 0].set_xlabel('Position along fiber')
    # axs[1, 1].set_xlabel('Position along fiber')
