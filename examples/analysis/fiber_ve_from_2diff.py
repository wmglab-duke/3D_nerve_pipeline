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
threshtarget = 'def'


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
    [5721],
    [5731],
    ['5Rdef'],
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
    for i, row in enumerate(threshdat.itertuples()):
        # plot 2D second differential
        inner = int(row.inner)
        fiber = int(row.fiber)
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

        axve.plot(v2, z_coords_loaded / 10000, alpha=0.1, color='b', label='2DEM' if i == 0 else '_')

        # run again for 3d
        inner3d = 0
        fiber3d = int(row.master_fiber_index)
        pve3 = os.path.join(base_n_sim3d, str(n_sim), 'data', 'inputs', f'inner{inner3d}_fiber{fiber3d}.dat')
        v3 = np.loadtxt(pve3)
        v3 = v3[1::11]

        # get z coordinate of each node
        zcoordpath = os.path.join(
            'samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'fibersets', str(n_sim), f'{fiber3d}.dat'
        )
        z_coords_arc_loaded = np.loadtxt(zcoordpath, skiprows=1)[:, 2][::11]

        fiberpath = os.path.join(
            'samples', str(samp3d), 'models', str(model), 'sims', '3', '3D_fiberset', f'{fiber3d}.dat'
        )
        fiber = nd_line(np.loadtxt(fiberpath, skiprows=1))
        loaded_3d = np.array([fiber.interp(d) for d in z_coords_arc_loaded])[:, 2]

        axve.plot(v3, loaded_3d / 10000, alpha=0.1, color='k', label='3DM' if i == 0 else '_')

    figve.set_size_inches(4, 6)
    leg = axve.legend()
    for lh in leg.legend_handles:
        lh.set_alpha(1)
    axve.set_ylabel('Distance along nerve (cm)')
    axve.set_xlabel('Electrical potential (V)')
    figve.savefig(f'test{nerve_label}.png')
