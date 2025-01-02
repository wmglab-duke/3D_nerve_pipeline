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
sim = 333
threshtarget = 'def'
early_exit = False
normalize = False

# List of n_sims to loop over
n_sims_list = [0, 5]  # Change as needed


def loadcoord(sim_object, sample, model, sim, n_sim, inner, fiber):
    """Load z-coordinates for the specified fiber."""
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


rho = [1, 1, 5]
for sample, samp3d, nerve_label in zip(
    [2521, 3721, 5721, 6721],
    [253, 373, 573, 653, 673],
    ['2Ldef', '3Rdef', '5Rdef', '6Rdef'],
):

    # Initialize figure and subplots for each nerve
    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey='row')
    # Axis references for clarity:
    #   axs[0, 0] -> Ve (2D)
    #   axs[0, 1] -> Ve (3D)
    #   axs[1, 0] -> 2nd Diff (2D)
    #   axs[1, 1] -> 2nd Diff (3D)

    sns.set_style('whitegrid')

    # For color mapping, one color per n_sim
    # colors = plt.cm.viridis(np.linspace(0, 1, len(n_sims_list)))
    colors = ['k', 'gray']

    for n_sim_i, color in zip(n_sims_list, colors):
        print(n_sim_i, sample)
        # Load the simulation object for this sample
        sim_object = Query.get_object(Object.SIMULATION, [sample, model, sim])

        # Load threshold data for this n_sim
        threshdat = pd.read_csv(f'thresh_unmatched_sim{sim}_{threshtarget}.csv')
        threshdat = threshdat.query(f'sample=={sample} & nsim=={n_sim_i} & model=={model} & sim=={sim}')
        threshdat.reset_index(drop=True, inplace=True)

        if len(threshdat) == 0:
            print(f'No data for sample={sample}, n_sim={n_sim_i}. Skipping.')
            continue

        # Base directories for 2D and 3D inputs
        base_n_sim = os.path.join('samples', str(sample), 'models', str(model), 'sims', str(sim), 'n_sims')
        base_n_sim3d = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'n_sims')

        for i, row in threshdat.iterrows():
            inner = int(row['inner'])
            fiber = int(row['fiber'])

            # -------- 2D Data --------
            pve2 = os.path.join(base_n_sim, str(n_sim_i), 'data', 'inputs', f'inner{inner}_fiber{fiber}.dat')
            v2 = np.loadtxt(pve2)
            # Remove first value and take every 11th value (myelinated fiber assumption)
            v2 = v2[1::11]

            # 2D second difference, normalized
            der2_2d = np.diff(np.diff(v2))
            if normalize:
                der2_2d /= np.max(np.abs(der2_2d))

            # Z-coords for 2D
            z_coords_loaded = loadcoord(sim_object, sample, model, sim, n_sim_i, inner, fiber)[::11]
            # For second difference, we have fewer points
            z_coords_der = (z_coords_loaded[1:] + z_coords_loaded[:-1]) / 2
            z_coords_der = (z_coords_der[1:] + z_coords_der[:-1]) / 2

            # Plot Ve (2D)
            axs[0, 0].plot(
                z_coords_loaded / 10000,  # distance on x-axis
                v2,  # Ve on y-axis
                # alpha=0.2,
                color=color,
                label=f'n_sim={n_sim_i}' if i == 0 else '',  # Label only once per subplot
            )

            # Plot 2nd Diff (2D)
            axs[1, 0].plot(
                z_coords_der / 10000,
                der2_2d,
                # alpha=0.2,
                color=color,
                label=f'n_sim={n_sim_i}' if i == 0 else '',
            )

            # -------- 3D Data --------
            fiber3d = int(row['master_fiber_index'])
            pve3 = os.path.join(base_n_sim3d, str(n_sim_i), 'data', 'inputs', f'inner0_fiber{fiber3d}.dat')
            v3 = np.loadtxt(pve3)
            v3 = v3[1::11]
            der3_3d = np.diff(np.diff(v3))
            if normalize:
                der3_3d /= np.max(np.abs(der3_3d))

            # Load 3D fiber coordinates
            zcoordpath = os.path.join(
                'samples',
                str(samp3d),
                'models',
                str(model),
                'sims',
                str(sim),
                'fibersets',
                str(n_sim_i),
                f'{fiber3d}.dat',
            )
            z_coords_arc_loaded = np.loadtxt(zcoordpath, skiprows=1)[:, 2][::11]
            z_coords_arc = (z_coords_arc_loaded[1:] + z_coords_arc_loaded[:-1]) / 2
            z_coords_arc = (z_coords_arc[1:] + z_coords_arc[:-1]) / 2

            # Use nd_line to convert arc distances to actual spatial coordinates
            fiberpath = os.path.join(
                'samples', str(samp3d), 'models', str(model), 'sims', '3', '3D_fiberset', f'{fiber3d}.dat'
            )
            fiber_obj = nd_line(np.loadtxt(fiberpath, skiprows=1))

            # Full node coords for Ve (3D) plot
            loaded_3d = np.array([fiber_obj.interp(d) for d in z_coords_arc_loaded])[:, 2]

            # For second difference plot
            z_coords3d = np.array([fiber_obj.interp(d) for d in z_coords_arc])[:, 2]

            # Plot Ve (3D)
            axs[0, 1].plot(
                loaded_3d / 10000,
                v3,
                # alpha=0.2,
                color=color,
                label=f'n_sim={n_sim_i}' if i == 0 else '',
            )

            # Plot 2nd Diff (3D)
            axs[1, 1].plot(
                z_coords3d / 10000,
                der3_3d,
                # alpha=0.2,
                color=color,
                label=f'n_sim={n_sim_i}' if i == 0 else '',
            )

            if early_exit and i > early_exit:
                break

    # ---- Final Labeling and Styling for this nerve_label ----
    axs[0, 0].set_title('2D Extrusion: $V_e$')
    axs[0, 1].set_title('3D: $V_e$')
    axs[1, 0].set_title('2D Extrusion: 2nd Diff')
    axs[1, 1].set_title('3D: 2nd Diff')

    # X labels
    axs[1, 0].set_xlabel('Distance (cm)')
    axs[1, 1].set_xlabel('Distance (cm)')

    # Y labels
    axs[0, 0].set_ylabel('Ve (V)')
    axs[1, 0].set_ylabel(f'Second Diff{" (normalized)" if normalize else ""}')

    # Add legends (one time per subplot)
    for ax_row in axs:
        for ax in ax_row:
            handles, labels = ax.get_legend_handles_labels()
            # Avoid repeating labels
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys(), loc='best')

    fig.suptitle(f'{nerve_label}: Ve & 2nd Diff (2D vs 3D) [n_sims={n_sims_list}]', fontsize=14)
    fig.tight_layout()

    outdir = 'plots/2diffnew'
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(
        os.path.join(outdir, f'{nerve_label}_{"-".join(map(str,n_sims_list))}_{sim}_2x2_norm{normalize}.png'), dpi=400
    )
    plt.close(fig)
