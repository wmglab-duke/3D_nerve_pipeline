import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

os.chdir('../../')
model = 0
sim = 3
nsim = 0
for sample, samp3d, nerve_label in zip(
    # [252, 372, 572, 652, 672], [253, 373, 573, 653, 673], ['2L', '3R', '5R', '6L', '6R']
    [252],
    [253],
    ['2L'],
):
    # todo: make loop sampples and save to disk
    threshdat = pd.read_csv('thresh_unmatched_sim3_og.csv')
    threshdat = threshdat.query(f'sample=={sample} & nsim==0 & model=={model} & sim=={sim}')
    threshdat.reset_index(drop=True, inplace=True)
    sns.set_style('whitegrid')
    # load in potentials for each fiber, calculate second differential, and plot
    base_n_sim = os.path.join('samples', str(sample), 'models', str(model), 'sims', str(sim), 'n_sims')
    base_n_sim3d = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'n_sims')
    # make figure with 2 subplots, share y axis
    fig, axs = plt.subplots(2, 1, sharey=True, sharex=True)

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
        # remove first value and take every 11th value for myelinated fibers
        v2 = v2[1::11]
        v3 = v3[1::11]
        # calculate second differential and normalize to maximum absolute value
        der2 = np.diff(np.diff(v2))
        der3 = np.diff(np.diff(v3))
        der2 = der2 / np.max(np.abs(der2))
        der3 = der3 / np.max(np.abs(der3))
        # reverse for 3d
        der3 = der3[::-1]
        # plot with transparency
        axs[0].plot(der2, alpha=0.1, color='r')
        axs[1].plot(der3, alpha=0.1, color='r')
        print(i)
    axs[0].set_ylabel('Normalized 2nd Differential')
    axs[1].set_ylabel('Normalized 2nd Differential')
    plt.xlabel('Position along fiber')
    plt.suptitle(f'2nd Differential of Ve for all {nerve_label} fibers')
    # label each axis
    axs[0].set_title('2D')
    axs[1].set_title('3D')
    plt.savefig(f'{nerve_label}_2diff.png')
