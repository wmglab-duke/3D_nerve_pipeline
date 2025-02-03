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

# List of n_sims to loop over
n_sims_list = [0, 5]  # Change as needed
mfc = 'none'

# 1) New option: single_fiber
#    - False => plot all fibers
#    - (fiber_number, nerve_label) => plot only that fiber within that nerve_label, with markers
single_fiber = False  # (10,'5Rdef')  # e.g. (2, '2Ldef') to restrict plotting


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


# Nerves and their 3D-sample counterparts
samples = [2521, 3721, 5721, 6721]
samples_3d = [2531, 3731, 5731, 6731]
nerve_labels = ['2Ldef', '3Rdef', '5Rdef', '6Rdef']

for sample, samp3d, nerve_label in zip(samples, samples_3d, nerve_labels):

    # If single_fiber is a tuple, only proceed if the nerve_label matches
    if single_fiber and single_fiber[1] != nerve_label:
        continue

    # --- Figure: 3 rows x 2 columns ---
    # Row 0 -> Ve (Extrusion vs True-3D)
    # Row 1 -> 2nd Diff (unnormalized)
    # Row 2 -> 2nd Diff (normalized)
    fig, axs = plt.subplots(3, 2, figsize=(10, 10), sharex=True, sharey='row')
    # Axis references for clarity:
    #  axs[0, 0] -> Ve (Extrusion)
    #  axs[0, 1] -> Ve (True-3D)
    #  axs[1, 0] -> 2nd Diff unnormalized (Extrusion)
    #  axs[1, 1] -> 2nd Diff unnormalized (True-3D)
    #  axs[2, 0] -> 2nd Diff normalized (Extrusion)
    #  axs[2, 1] -> 2nd Diff normalized (True-3D)

    sns.set(style='white', context='paper')

    # Predefine colors for each n_sim
    colors = ['gray', 'k']  # or use a colormap if more n_sims

    for n_sim_i, color in zip(n_sims_list, colors):
        # Create a more informative label:
        if n_sim_i == 0:
            sim_label = '3 μm fibers'
        elif n_sim_i == 5:
            sim_label = '13 μm fibers'
        else:
            sim_label = f'n_sim={n_sim_i}'

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

        # Base directories for 2D (extrusion) and 3D inputs
        base_n_sim = os.path.join('samples', str(sample), 'models', str(model), 'sims', str(sim), 'n_sims')
        base_n_sim3d = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), 'n_sims')

        # Loop over each fiber row in the threshold data
        for i, row in threshdat.iterrows():
            inner = int(row['inner'])
            fiber = int(row['fiber'])

            # If single_fiber is set, skip unless it matches the (fiber_number, nerve_label) we want
            if single_fiber:
                if i != single_fiber[0]:
                    continue
                # We'll add markers in that case
                marker_style = 'o'
            else:
                marker_style = None

            # -------- Extrusion (2D) Data --------
            pve2_path = os.path.join(base_n_sim, str(n_sim_i), 'data', 'inputs', f'inner{inner}_fiber{fiber}.dat')
            v2 = np.loadtxt(pve2_path)
            # Remove first value and take every 11th value (myelinated fiber assumption)
            v2 = v2[1::11]

            # 2D second difference (unnormalized and normalized)
            der2_2d_unnorm = np.diff(np.diff(v2))
            max_der2_2d = np.max(np.abs(der2_2d_unnorm)) if len(der2_2d_unnorm) > 0 else 1
            der2_2d_norm = der2_2d_unnorm / max_der2_2d if max_der2_2d != 0 else der2_2d_unnorm

            # Z-coords for 2D
            z_coords_loaded = loadcoord(sim_object, sample, model, sim, n_sim_i, inner, fiber)[::11]

            # For second difference, we have fewer points
            # (the midpoint of adjacent points for the first diff, then midpoint again for second diff)
            z_coords_der = (z_coords_loaded[1:] + z_coords_loaded[:-1]) / 2
            z_coords_der = (z_coords_der[1:] + z_coords_der[:-1]) / 2

            # -------- Plot Extrusion (2D) Subplots --------
            axs[0, 0].plot(
                z_coords_loaded / 10000,  # distance on x-axis
                v2,  # Ve on y-axis
                color=color,
                marker=marker_style,
                markerfacecolor=mfc,
                linestyle='-',
                label=sim_label if i == 0 else '',
            )
            axs[1, 0].plot(
                z_coords_der / 10000,
                der2_2d_unnorm,
                color=color,
                marker=marker_style,
                markerfacecolor=mfc,
                linestyle='-',
                label=sim_label if i == 0 else '',
            )
            axs[2, 0].plot(
                z_coords_der / 10000,
                der2_2d_norm,
                color=color,
                marker=marker_style,
                markerfacecolor=mfc,
                linestyle='-',
                label=sim_label if i == 0 else '',
            )

            # -------- True-3D Data --------
            fiber3d = int(row['master_fiber_index'])  # from your CSV data
            pve3_path = os.path.join(base_n_sim3d, str(n_sim_i), 'data', 'inputs', f'inner0_fiber{fiber3d}.dat')
            v3 = np.loadtxt(pve3_path)
            v3 = v3[1::11]

            der3_3d_unnorm = np.diff(np.diff(v3))
            max_der3_3d = np.max(np.abs(der3_3d_unnorm)) if len(der3_3d_unnorm) > 0 else 1
            der3_3d_norm = der3_3d_unnorm / max_der3_3d if max_der3_3d != 0 else der3_3d_unnorm

            # Load 3D fiber coordinates for arc length
            zcoordpath_3d = os.path.join(
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
            z_coords_arc_loaded = np.loadtxt(zcoordpath_3d, skiprows=1)[:, 2][::11]
            z_coords_arc = (z_coords_arc_loaded[1:] + z_coords_arc_loaded[:-1]) / 2
            z_coords_arc = (z_coords_arc[1:] + z_coords_arc[:-1]) / 2

            # Use nd_line to convert arc distances to spatial coords (particularly the z-component in this example)
            fiberpath_3d = os.path.join(
                'samples', str(samp3d), 'models', str(model), 'sims', '3', '3D_fiberset', f'{fiber3d}.dat'
            )
            fiber_obj_3d = nd_line(np.loadtxt(fiberpath_3d, skiprows=1))

            # Full node coords for Ve (3D) plot
            loaded_3d_z = np.array([fiber_obj_3d.interp(d) for d in z_coords_arc_loaded])[:, 2]
            # 2nd diff subplot coords
            z_coords3d = np.array([fiber_obj_3d.interp(d) for d in z_coords_arc])[:, 2]

            # -------- Plot True-3D Subplots --------
            axs[0, 1].plot(
                loaded_3d_z / 10000,
                v3,
                color=color,
                marker=marker_style,
                markerfacecolor=mfc,
                linestyle='-',
                label=sim_label if i == 0 else '',
            )
            axs[1, 1].plot(
                z_coords3d / 10000,
                der3_3d_unnorm,
                color=color,
                marker=marker_style,
                markerfacecolor=mfc,
                linestyle='-',
                label=sim_label if i == 0 else '',
            )
            axs[2, 1].plot(
                z_coords3d / 10000,
                der3_3d_norm,
                color=color,
                marker=marker_style,
                markerfacecolor=mfc,
                linestyle='-',
                label=sim_label if i == 0 else '',
            )

            if early_exit and i > early_exit:
                break

    # ---- Final Labeling and Styling for this nerve_label ----
    axs[0, 0].set_title('Extrusion')
    axs[0, 1].set_title('True-3D')

    # Legend
    # create handles
    handles = [
        plt.Line2D([0], [0], color='gray', marker=marker_style, linestyle='-', markersize=5, label='3', mfc=mfc),
        plt.Line2D([0], [0], color='k', marker=marker_style, linestyle='-', markersize=5, label='13', mfc=mfc),
    ]
    axs[0, 0].legend(handles=handles, title='D: (μm)', loc='upper right')

    # X labels (bottom row only)
    axs[2, 0].set_xlabel('Distance (cm)')
    axs[2, 1].set_xlabel('Distance (cm)')

    # Y labels
    axs[0, 0].set_ylabel('Ve (V)')
    axs[1, 0].set_ylabel('2nd Diff.')
    axs[2, 0].set_ylabel('Normalized 2nd Diff.')

    if not single_fiber:
        fig.suptitle(f'{nerve_label}: Ve & 2nd Diff (Extrusion vs. True-3D)', fontsize=14)
    plt.gcf().set_size_inches(6, 6)
    fig.tight_layout()

    outdir = 'plots/2diffnew'
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, f'{nerve_label}_{"-".join(map(str,n_sims_list))}_{sim}_3rows.png'), dpi=400)
    plt.show()
    plt.close(fig)
