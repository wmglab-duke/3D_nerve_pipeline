import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

os.chdir('../../')
sample = 252
samp3d = 253
model = 0
sim = 3
threshdat = pd.read_csv('thresh.csv')
threshdat = threshdat.query(f'sample=={sample} & nsim==0 & model=={model} & sim=={sim}')
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
    pve2 = os.path.join(base_n_sim, str(0), 'data', 'inputs', f'inner{inner}_fiber{fiber}.dat')
    pve3 = os.path.join(base_n_sim3d, str(0), 'data', 'inputs', f'inner{inner3d}_fiber{fiber3d}.dat')
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
axs[0].set_ylabel('Normalized 2nd Differential')
axs[1].set_ylabel('Normalized 2nd Differential')
plt.xlabel('Position along fiber')
plt.suptitle('2nd Differential of Ve for all 2L fibers')
# label each axis
axs[0].set_title('2D')
axs[1].set_title('3D')
plt.show()
