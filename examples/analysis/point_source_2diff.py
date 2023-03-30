import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from nd_line.nd_line import nd_line

os.chdir('../../')
sample = 252
samp3d = 253
nerve_label = '2L'
model = 0
sim = 3

def pt_src(sourcepoint, nodepoint,current = 1,rho=1):
    # sourcepoint = (x,y,z) of source
    # nodepoint = (x,y,z) of node
    # current = current in source
    # rho = resistivity of medium
    # returns potential at node

    # calculate distance between source and node
    dist = np.sqrt((sourcepoint[0]-nodepoint[0])**2 + (sourcepoint[1]-nodepoint[1])**2 + (sourcepoint[2]-nodepoint[2])**2)
    # calculate potential
    v = current/(4*np.pi*rho*dist)
    return v

for sample, samp3d, nerve_label in zip(
    [252, 372, 572, 652, 672], [253, 373, 573, 653, 673], ['2L', '3R', '5R', '6L', '6R']
):
    # todo: make loop samples and save to disk
    threshdat = pd.read_csv('thresh_unmatched_sim3.csv')
    threshdat = threshdat.query(f'sample=={sample} & nsim==0 & model=={model} & sim=={sim}')
    threshdat.reset_index(drop=True, inplace=True)
    sns.set_style('whitegrid')
    # load in potentials for each fiber, calculate second differential, and plot
    base_n_sim = os.path.join('samples', str(sample), 'models', str(model), 'sims', str(sim), 'n_sims')
    base_n_sim3d = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'n_sims')
    # make figure with 2 subplots, share y axis
    fig, axs = plt.subplots(3, 1, sharey=True, sharex=True)

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
        # print(i)
        #now load the 3D fiber and calculate the potential at each node
        fiberpath = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), '3D_fiberset', f'{fiber3d}.dat')
        fiber = nd_line(np.loadtxt(fiberpath,skiprows = 1))
        #assert that z coords are monotonic
        try:assert np.all(np.diff(fiber.points[:,2])<0)
        except AssertionError:print(i,'is not monotonic')
        #load in fiber coordinates from ascent fiber
        fiber3d = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'fibersets', str(0),f'{fiber3d}.dat')
        dcoords = np.loadtxt(fiber3d,skiprows=1)[:,2][1::11] #arc distances along line
        #calculate potential at each node
        v = np.zeros(len(dcoords))
        source_point1 = (-6.394884621840902E-13,
1696.9375429688143,
21050.0
)
        source_point2 = (2.2737367544323206E-13,
1430.8222093677691,
29050.000000000004

)
        for j, d in enumerate(dcoords):
            node_point = fiber.interp(d)
            v[j] = pt_src(source_point1,node_point)
            v[j]+=pt_src(source_point2,node_point,current=-1)
        #calculate second differential and normalize to maximum absolute value
        der2 = np.diff(np.diff(v))
        der2 = der2 / np.max(np.abs(der2))
        #plot with transparency
        axs[2].plot(der2, alpha=0.1, color='r')
    axs[1].set_ylabel('Normalized 2nd Differential')
    plt.xlabel('Node Number')
    plt.suptitle(f'2nd Differential of Ve for all {nerve_label} fibers')
    # label each axis
    axs[0].set_title('2D')
    axs[1].set_title('3D')
    axs[2].set_title('3D (point source, homogeneous isotropic medium)')
    fig.set_size_inches(6,6)
    plt.savefig(f'{nerve_label}_2diff.png')
