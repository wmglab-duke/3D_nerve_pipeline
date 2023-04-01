import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from nd_line.nd_line import nd_line

os.chdir('../../')

from src.core import Query
from src.utils import Object

sample = 252
samp3d = 253
nerve_label = '2L'
model = 0
sim = 3


def pt_src(sourcepoint, nodepoint, current=1, rho=[1,1,1]):
    # sourcepoint = (x,y,z) of source
    # nodepoint = (x,y,z) of node
    # current = current in source
    # rho = resistivity of medium
    # returns potential at node

    # calculate distance between source and node
    dists = np.array(sourcepoint) - np.array(nodepoint)
    # calculate potential
    v = current / (4 * np.pi * np.sqrt(rho[1]*rho[2]*dists[0]**2 + rho[0]*rho[2]*dists[1]**2 + rho[0]*rho[1]*dists[2]**2))
    return v

def loadcoord(sim_object,sample,model,sim,n_sim):
    # load the corresponding fiber coordinates
    for t, (p_i, _) in enumerate(sim_object.master_product_indices):
        if t == n_sim:
            potentials_ind = p_i
            break

    active_src_ind, fiberset_ind = sim_object.potentials_product[potentials_ind]
    master_fiber_ind = sim_object.indices_n_to_fib(fiberset_index=fiberset_ind, inner_index=inner,
                                                   local_fiber_index=fiber)

    fiber_coords_path = os.path.join(
        'samples',
        str(sample),
        'models',
        str(model),
        'sims',
        str(sim),
        'fibersets',
        str(fiberset_ind),
        f'{master_fiber_ind}.dat',
    )
    z_coords = np.loadtxt(fiber_coords_path, skiprows=1)[:, 2]
    return z_coords

n_sim=0
rho = [1,1,10]
for sample, samp3d, nerve_label in zip(
    [252], [253], ['2L']
):
    sim_object = Query.get_object(Object.SIMULATION, [sample, model, sim])

    # todo: make loop samples and save to disk
    threshdat = pd.read_csv('thresh_unmatched_sim3.csv')
    threshdat = threshdat.query(f'sample=={sample} & nsim==0 & model=={model} & sim=={sim}')
    threshdat.reset_index(drop=True, inplace=True)
    sns.set_style('whitegrid')
    # load in potentials for each fiber, calculate second differential, and plot
    base_n_sim = os.path.join('samples', str(sample), 'models', str(model), 'sims', str(sim), 'n_sims')
    base_n_sim3d = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'n_sims')
    # make figure with 2 subplots, share y axis
    fig, axs = plt.subplots(1, 3, sharey=True, sharex=True)

    # loop through each fiber
    for i, row in threshdat.iterrows():
        #plot 2D second differential
        inner = int(row['inner'])
        fiber = int(row['fiber'])
        pve2 = os.path.join(base_n_sim, str(n_sim), 'data', 'inputs', f'inner{inner}_fiber{fiber}.dat')
        v2 = np.loadtxt(pve2)
        print(i)

        # remove first value and take every 11th value for myelinated fibers
        v2 = v2[1::11]

        # calculate second differential and normalize to maximum absolute value
        der2 = np.diff(np.diff(v2))
        der2 = der2 / np.max(np.abs(der2))

        #get z coordinate of each node
        z_coords = loadcoord(sim_object,sample,model,sim,n_sim)[::11][2:]
        zspacing = np.diff(z_coords)[0]

        # plot with transparency
        axs[0].plot(der2,z_coords, alpha=0.1, color='r')

        #run again for 3d
        inner3d = 0
        fiber3d = int(row['master_fiber_index'])
        pve3 = os.path.join(base_n_sim3d, str(n_sim), 'data', 'inputs', f'inner{inner3d}_fiber{fiber3d}.dat')
        v3 = np.loadtxt(pve3)
        v3 = v3[1::11]
        der3 = np.diff(np.diff(v3))
        der3 = der3 / np.max(np.abs(der3))
        # reverse for 3d
        # der3 = der3[::-1]

        #get z coordinate of each node
        zcoordpath = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'fibersets',str(n_sim), f'{fiber3d}.dat')
        z_coords_arc = np.loadtxt(zcoordpath, skiprows=1)[:,2][::11][2:]
        fiberpath = os.path.join(
            'samples', str(samp3d), 'models', str(model), 'sims', str(sim), '3D_fiberset', f'{fiber3d}.dat'
        )
        fiber = nd_line(np.loadtxt(fiberpath, skiprows=1))
        z_coords3d = np.array([fiber.interp(d) for d in z_coords_arc])[:,2]

        axs[1].plot(der3, z_coords3d,alpha=0.1, color='r')

        #plot 3D second differential for supersamples
        # load in fiber coordinates from ascent fiber
        # sspath0 = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'ss_bases','0', f'{fiber3d}.dat')
        # sspath1 = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'ss_bases','1', f'{fiber3d}.dat')

        # ss1 = np.loadtxt(sspath0, skiprows=1)
        # ss2 = np.loadtxt(sspath1, skiprows=1)

        # superve = -ss1+ss2
        
        # superpoints = fiber.points[::100][:,2]
        
        # superdiff=np.diff(np.diff(superve[::100]))

        # plt.plot(superpoints[2:],superdiff, alpha=0.1, color='r')


        # now load the 3D fiber and calculate the point source potential at each node
        # assert that z coords are monotonic
        try:
            assert np.all(np.diff(fiber.points[:, 2]) < 0)
        except AssertionError:
            print(i, 'is not monotonic')
        # load in fiber coordinates from ascent fiber
        fiber3dpath = os.path.join(
            'samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'fibersets', str(n_sim), f'{fiber3d}.dat'
        )
        dcoords = np.loadtxt(fiber3dpath, skiprows=1)[:,2][1::11]  # arc distances along line
        # calculate potential at each node
        v = np.zeros(len(dcoords))
        source_point1 = (-6.394884621840902e-13, 1696.9375429688143, 21050.0)
        source_point2 = (2.2737367544323206e-13, 1430.8222093677691, 29050.000000000004)
        zs = []
        for j, d in enumerate(dcoords):
            node_point = fiber.interp(d)
            v[j] = pt_src(source_point1, node_point, rho=rho,current=-1)
            v[j] += pt_src(source_point2, node_point,rho=rho)
            zs.append(node_point[2])
        # calculate second differential and normalize to maximum absolute value
        der2 = np.diff(np.diff(v))
        der2 = der2 / np.max(np.abs(der2))
        # plot with transparency
        axs[2].plot(der2, zs[2:],alpha=0.1, color='r')
        #find maximum z coordinate
    #verical lines to each axis at z coords of point sources
    axs[0].axhline(y=source_point1[2], color='k', linestyle='--')
    axs[0].axhline(y=source_point2[2], color='k', linestyle='--')
    axs[1].axhline(y=source_point1[2], color='k', linestyle='--')
    axs[1].axhline(y=source_point2[2], color='k', linestyle='--')
    axs[2].axhline(y=source_point1[2], color='k', linestyle='--')
    axs[2].axhline(y=source_point2[2], color='k', linestyle='--',label='source position')
    zmax = np.max(fiber.points[:, 2])
    # find and plot voltage along a fiber at [0,0,0] to [0,0,zmax]
    for j, z in enumerate(np.arange(0, zmax, zspacing)):
        node_point = np.array([0, 0, z])
        v[j] = pt_src(source_point1, node_point,rho=rho, current=-1)
        v[j] += pt_src(source_point2, node_point,rho= rho)
    # calculate second differential and normalize to maximum absolute value
    der2 = np.diff(np.diff(v))
    der2 = der2 / np.max(np.abs(der2))
    # plot with transparency
    axs[2].plot(der2,np.arange(0, zmax, zspacing)[2:], alpha=1, color='b',linestyle='--',label='2D-path')
    # label axes
    axs[0].set_ylabel('z-coordinate (μm)')
    axs[1].set_xlabel('Normalized Second Differential')
    plt.suptitle(f'2nd Differential of Ve for all {nerve_label} fibers')
    # label each axis
    axs[0].set_title('2D')
    axs[1].set_title('3D')
    axs[2].set_title('3D (point source, \nhomogeneous \nanisotropic medium)')
    plt.subplots_adjust(top=0.8)
    plt.legend()
    plt.savefig(f'{nerve_label}_2diff.png',dpi=400)
