import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from nd_line.nd_line import nd_line

os.chdir('../../')

from src.core import Query
from src.utils import Object

model = 0
sim = 3
n_sim = 0
threshtarget = 'og'


def pt_src(sourcepoint, nodepoint, current=1, rho=[1, 1, 1]):
    # sourcepoint = (x,y,z) of source
    # nodepoint = (x,y,z) of node
    # current = current in source
    # rho = resistivity of medium
    # returns potential at node

    # calculate distance between source and node
    dists = np.array(sourcepoint) - np.array(nodepoint)
    # calculate potential
    v = current / (
        4
        * np.pi
        * np.sqrt(rho[1] * rho[2] * dists[0] ** 2 + rho[0] * rho[2] * dists[1] ** 2 + rho[0] * rho[1] * dists[2] ** 2)
    )
    return v


def loadcoord(sim_object, sample, model, sim, n_sim):
    # load the corresponding fiber coordinates
    for t, (p_i, _) in enumerate(sim_object.master_product_indices):
        if t == n_sim:
            potentials_ind = p_i
            break

    active_src_ind, fiberset_ind = sim_object.potentials_product[potentials_ind]
    master_fiber_ind = sim_object.indices_n_to_fib(
        fiberset_index=fiberset_ind, inner_index=inner, local_fiber_index=fiber
    )

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


# TODO get peak 2diffs for each sample and compare their location and value. start with 2L, will need to get contact centers for the rest of them
rho = [1, 1, 5]
for sample, samp3d, nerve_label in zip(
    # [252, 372, 572, 652, 672], [253, 373, 573, 653, 673], ['2L', '3R', '5R', '6L', '6R']
    # [252],
    # [253],
    # ['2L'],
    [2521, 3721, 5721, 6721],
    [253, 373, 573, 653, 673],
    ['2Ldef', '3Rdef', '5Rdef', '6Rdef'],
    # [2529, 3729, 5729, 6729], [2531, 3731, 5731, 6531, 6731], ['2Lasc', '3Rasc', '5Rasc', '6Rasc']
):
    sim_object = Query.get_object(Object.SIMULATION, [sample, model, sim])

    # todo: make loop samples and save to disk
    threshdat = pd.read_csv(f'thresh_unmatched_sim{sim}_{threshtarget}.csv')
    threshdat = threshdat.query(f'sample=={sample} & nsim=={n_sim} & model=={model} & sim=={sim}')
    threshdat.reset_index(drop=True, inplace=True)
    sns.set_style('whitegrid')
    # load in potentials for each fiber, calculate second differential, and plot
    base_n_sim = os.path.join('samples', str(sample), 'models', str(model), 'sims', str(sim), 'n_sims')
    base_n_sim3d = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'n_sims')
    # make figure with 2 subplots, share y axis
    fig, axs = plt.subplots(1, 3, sharey=True, sharex=True)

    figve, axve = plt.subplots()

    assert len(threshdat) > 0

    # loop through each fiber
    for i, row in threshdat.iterrows():
        # plot 2D second differential
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

        # get z coordinate of each node
        z_coords_loaded = loadcoord(sim_object, sample, model, sim, n_sim)
        z_coords_loaded = z_coords_loaded[::11]
        zspacing = np.diff(z_coords_loaded)[0]
        # get new point for second difference
        z_coords = z_coords_loaded.copy()
        z_coords = (z_coords[1:] + z_coords[:-1]) / 2
        z_coords = (z_coords[1:] + z_coords[:-1]) / 2

        # plot with transparency
        axs[0].plot(der2, z_coords / 10000, alpha=0.1, color='r')

        axve.plot(v2, z_coords_loaded / 10000, alpha=0.1, color='b', label='2DEM' if i == 0 else '_')

        # run again for 3d
        inner3d = 0
        fiber3d = int(row['master_fiber_index'])
        pve3 = os.path.join(base_n_sim3d, str(n_sim), 'data', 'inputs', f'inner{inner3d}_fiber{fiber3d}.dat')
        v3 = np.loadtxt(pve3)
        v3 = v3[1::11]
        der3 = np.diff(np.diff(v3))
        der3 = der3 / np.max(np.abs(der3))
        # reverse for 3d
        # der3 = der3[::-1]

        # get z coordinate of each node
        zcoordpath = os.path.join(
            'samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'fibersets', str(n_sim), f'{fiber3d}.dat'
        )
        z_coords_arc_loaded = np.loadtxt(zcoordpath, skiprows=1)[:, 2][::11]
        z_coords_arc = z_coords_arc_loaded.copy()
        # get new points for second difference
        z_coords_arc = (z_coords_arc[1:] + z_coords_arc[:-1]) / 2
        z_coords_arc = (z_coords_arc[1:] + z_coords_arc[:-1]) / 2

        fiberpath = os.path.join(
            'samples', str(samp3d), 'models', str(model), 'sims', '3', '3D_fiberset', f'{fiber3d}.dat'
        )
        fiber = nd_line(np.loadtxt(fiberpath, skiprows=1))
        z_coords3d = np.array([fiber.interp(d) for d in z_coords_arc])[:, 2]
        loaded_3d = np.array([fiber.interp(d) for d in z_coords_arc_loaded])[:, 2]

        axs[1].plot(der3, z_coords3d / 10000, alpha=0.1, color='r')
        axve.plot(v3, loaded_3d / 10000, alpha=0.1, color='k', label='3DM' if i == 0 else '_')

        try:
            assert np.all(np.diff(fiber.points[:, 2]) < 0)
        except AssertionError:
            print(i, 'is not monotonic')
        # load in fiber coordinates from ascent fiber
        fiber3dpath = os.path.join(
            'samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'fibersets', str(n_sim), f'{fiber3d}.dat'
        )
        dcoords = np.loadtxt(fiber3dpath, skiprows=1)[:, 2][1::11]  # arc distances along line
        # calculate potential at each node
        v = np.zeros(len(dcoords))
        # load source points
        pcsdir = os.path.join('input', 'pcs_locs')
        source_point1 = np.loadtxt(os.path.join(pcsdir, f'{nerve_label}DS5', 'pcs1.dat'))
        source_point2 = np.loadtxt(os.path.join(pcsdir, f'{nerve_label}DS5', 'pcs2.dat'))
        zs = []
        for j, d in enumerate(dcoords):
            node_point = fiber.interp(d)
            v[j] = pt_src(source_point1, node_point, rho=rho, current=-1)
            v[j] += pt_src(source_point2, node_point, rho=rho)
            zs.append(node_point[2])
        zs = np.array(zs)
        zs = (zs[1:] + zs[:-1]) / 2
        zs = (zs[1:] + zs[:-1]) / 2
        # calculate second differential and normalize to maximum absolute value
        der2 = np.diff(np.diff(v))
        der2 = der2 / np.max(np.abs(der2))
        # plot with transparency
        axs[2].plot(der2, np.array(zs) / 10000, alpha=0.1, color='r')
        # find maximum z coordinate
    # verical lines to each axis at z coords of point sources
    axs[0].axhline(y=source_point1[2] / 10000, color='k', linestyle='--')
    axs[0].axhline(y=source_point2[2] / 10000, color='k', linestyle='--')
    axs[1].axhline(y=source_point1[2] / 10000, color='k', linestyle='--')
    axs[1].axhline(y=source_point2[2] / 10000, color='k', linestyle='--')
    axs[2].axhline(y=source_point1[2] / 10000, color='k', linestyle='--')
    axs[2].axhline(y=source_point2[2] / 10000, color='k', linestyle='--', label='source position')
    zmax = np.max(fiber.points[:, 2])
    # find and plot voltage along a fiber at [0,0,0] to [0,0,zmax]
    v2 = np.zeros(len(np.arange(0, zmax, zspacing)))
    for j, z in enumerate(np.arange(0, zmax, zspacing)):
        node_point = np.array([0, 0, z])
        v2[j] = pt_src(source_point1, node_point, rho=rho, current=-1)
        v2[j] += pt_src(source_point2, node_point, rho=rho)
    # calculate second differential and normalize to maximum absolute value
    der2 = np.diff(np.diff(v2))
    der2 = der2 / np.max(np.abs(der2))
    # plot with transparency
    thisz = np.arange(0, zmax, zspacing)
    thisz = (thisz[1:] + thisz[:-1]) / 2
    thisz = (thisz[1:] + thisz[:-1]) / 2

    axs[2].plot(der2, thisz / 10000, alpha=1, color='b', linestyle='--', label='2D-path')
    # label axes
    axs[0].set_ylabel('z-coordinate (cm)')
    axs[1].set_xlabel('Normalized Second Differential')
    fig.suptitle(f'2nd Differential of Ve for all {nerve_label} fibers')
    # label each axis
    axs[0].set_title('2D')
    axs[1].set_title('3D')
    axs[2].set_title('3D (point sources)')
    fig.subplots_adjust(top=0.8)
    fig.axes[-1].legend()
    fig.savefig(f'plots/2diff/{nerve_label}_{n_sim}_{sim}2diff.png', dpi=400)
    figve.set_size_inches(4, 6)
    axve.legend()
    axve.set_ylabel('Distance along nerve (cm)')
    axve.set_xlabel('Electrical potential (V)')
    figve.savefig(f'plots/2diff/{nerve_label}_{n_sim}_{sim}ve.png', dpi=400)
