import os
import pickle
import random
import sys

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from shapely.geometry import Point
from shapely.strtree import STRtree

font = {'size': 12}
matplotlib.rc('font', **font)

# If needed, adjust your NEURON path
sys.path.append(r'C:\nrn\lib\python')  # noqa: E800

# Move up two directories so script can find relative paths (adjust if needed)
os.chdir('../..')

# Palettes for 2D vs 3D
pal2d3d = ['#d95f02', '#7570b3']

# Dictionary of sample IDs for undeformed (ud) and deformed (def)
# Fill in the correct IDs for '2L','3R','6R' as needed
sample_ids = {
    '2L': {
        '2d_ud': None,  # e.g. 211
        '3d_ud': None,  # e.g. 213
        '2d_def': None,  # e.g. 2111
        '3d_def': None,  # e.g. 2131
    },
    '3R': {
        '2d_ud': 372,
        '3d_ud': 373,
        '2d_def': 3721,
        '3d_def': 3731,
    },
    '5R': {
        '2d_ud': 572,
        '3d_ud': 573,
        '2d_def': 5721,
        '3d_def': 5731,
    },
    '6R': {
        '2d_ud': None,
        '3d_ud': None,
        '2d_def': None,
        '3d_def': None,
    },
}

# We will sample at multiple z-step intervals
z_steps = [500, 1000]  # adjust as you like


# Helper function to load + compute for given sample, 2D/3D, and deformation
def load_fiber_data(sample_name: str, samp2d: int, samp3d: int):
    """Loads data for a single sample (2D & 3D for undeformed or deformed),
    computes Ve and second difference, and returns a dictionary."""
    # Directory containing fiber files
    fiberpath = os.path.join(os.getcwd(), fr'samples\{samp3d}\models\0\sims\3\3D_fiberset')

    # Load all fibers from .dat files
    fibers = {}
    for file in os.listdir(fiberpath):
        if file.endswith('.dat'):
            fibers[int(file.replace('.dat', ''))] = np.loadtxt(os.path.join(fiberpath, file), skiprows=1)

    data = {
        key: {
            "z_coords": arr[:, 2],
            "3d_coords": arr,
            "ve": None,
            "ve_2diff": None,
            "z_coords_2d": None,
            "ve_2d": None,
            "ve_2diff_2d": None,
        }
        for key, arr in fibers.items()
    }

    # For each fiber, load the Ve from both contacts (for 3D)
    for mfi, value in data.items():
        contacts = [0, 1]
        weights = [1, -1]

        # ============ True 3D ============
        ves_3d = []
        for contact, weight in zip(contacts, weights):
            threedpath = rf'D:\threed_final\samples\{samp3d}\models\0\sims\3\ss_bases\{contact}\{mfi}.dat'
            # multiply by 1000 for mV, multiply by weight to handle polarity
            ves_3d.append(np.loadtxt(threedpath, skiprows=1) * 1000 * weight)
        value['ve'] = np.sum(ves_3d, axis=0)
        # second derivative
        value['ve_2diff'] = np.diff(np.diff(value['ve']))

        # ============ 2D Extrusion ============
        ves_2d = []
        for contact, weight in zip(contacts, weights):
            twoodpath = rf'D:\threed_final\samples\{samp2d}\models\0\sims\3\ss_bases\{contact}\{mfi}.dat'
            ves_2d.append(np.loadtxt(twoodpath, skiprows=1) * 1000 * weight)
        value['ve_2d'] = np.sum(ves_2d, axis=0)
        value['ve_2diff_2d'] = np.diff(np.diff(value['ve_2d']))

        # 2D z-coords
        twoodc = rf'D:\threed_final\samples\{samp2d}\models\0\sims\3\ss_coords\{mfi}.dat'
        value['z_coords_2d'] = np.loadtxt(twoodc, skiprows=1)[:, 2]

    return data


def get_alldata_for_zsteps(
    sample_name: str, samp2d_ud: int, samp3d_ud: int, samp2d_def: int, samp3d_def: int, zstep: int
):
    """Loads data for both undeformed and deformed sets for a single zstep,
    returns (df_Ve, df_2diff) after sampling along z."""
    # Load undeformed
    data_og = load_fiber_data(sample_name, samp2d_ud, samp3d_ud)
    # Load deformed
    data_def = load_fiber_data(sample_name, samp2d_def, samp3d_def)

    # We'll sample from z = 50 to 50,000 in increments of zstep
    zs = np.arange(50, 50100, zstep)

    # Build the data frames for Ve and second difference
    alldata = {}
    all_2diff_data = {}

    for label, dataset in zip(['og', 'def'], [data_og, data_def]):
        for mfi, value in dataset.items():
            # sample the 2D/3D Ve at the chosen zs
            # 2D
            ve_2_inds = [np.argmin(np.abs(value['z_coords_2d'] - z)) for z in zs]
            ve_2d_at_zs = value['ve_2d'][ve_2_inds]
            # 3D
            ve_3_inds = [np.argmin(np.abs(value['z_coords'] - z)) for z in zs]
            ve_3d_at_zs = value['ve'][ve_3_inds]

            # Save back
            dataset[mfi]['ve_2d_at_zs'] = ve_2d_at_zs
            dataset[mfi]['ve_3d_at_zs'] = ve_3d_at_zs

            # Second difference at those sample points
            # Because these are new samples, let's do np.diff on them
            dataset[mfi]['ve_2diff_2d_at_zs'] = np.diff(np.diff(ve_2d_at_zs))
            dataset[mfi]['ve_2diff_3d_at_zs'] = np.diff(np.diff(ve_3d_at_zs))

        # Convert to DataFrame for Ve
        df_2d = pd.DataFrame({key: val['ve_2d_at_zs'] for key, val in dataset.items()})
        df_3d = pd.DataFrame({key: val['ve_3d_at_zs'] for key, val in dataset.items()})
        df_2d.index = zs
        df_3d.index = zs
        df_2d = df_2d.reset_index().melt(id_vars='index')
        df_3d = df_3d.reset_index().melt(id_vars='index')
        df_2d.columns = ['z', 'fiber', 've']
        df_3d.columns = ['z', 'fiber', 've']
        df_2d['dim'] = '2d'
        df_3d['dim'] = '3d'
        df = pd.concat([df_2d, df_3d])
        df['deformed'] = 'no' if label == 'og' else 'yes'

        # Save in dict
        if label not in alldata:
            alldata[label] = df
        else:
            # or just store them in a list to combine later
            pass

        # Convert to DataFrame for 2nd diff
        df_2d_2diff = pd.DataFrame({key: val['ve_2diff_2d_at_zs'] for key, val in dataset.items()})
        df_3d_2diff = pd.DataFrame({key: val['ve_2diff_3d_at_zs'] for key, val in dataset.items()})

        # Because we did np.diff(np.diff(...)), the new index is one smaller by 2
        # i.e. from zs[1:-1]
        if len(zs) > 2:
            df_2d_2diff.index = zs[1:-1]
            df_3d_2diff.index = zs[1:-1]
        else:
            # If your step array is very small, watch out for edge cases
            df_2d_2diff.index = []
            df_3d_2diff.index = []

        df_2d_2diff = df_2d_2diff.reset_index().melt(id_vars='index')
        df_3d_2diff = df_3d_2diff.reset_index().melt(id_vars='index')
        df_2d_2diff.columns = ['z', 'fiber', 've']
        df_3d_2diff.columns = ['z', 'fiber', 've']
        df_2d_2diff['dim'] = '2d'
        df_3d_2diff['dim'] = '3d'
        df2 = pd.concat([df_2d_2diff, df_3d_2diff])
        df2['deformed'] = 'no' if label == 'og' else 'yes'

        if label not in all_2diff_data:
            all_2diff_data[label] = df2

    # Combine undeformed + deformed for Ve
    combined_alldata = pd.concat(list(alldata.values()), ignore_index=True)
    # Combine undeformed + deformed for 2diff
    combined_2diff_data = pd.concat(list(all_2diff_data.values()), ignore_index=True)

    return combined_alldata, combined_2diff_data


# ==============================================================
# Main script: loop over samples, z-steps, then plot subplots
# ==============================================================

samples = ['2L', '3R', '5R', '6R']

for samp in samples:
    # Unpack the correct IDs for the sample
    ids = sample_ids[samp]
    if any(v is None for v in ids.values()):
        print(f"Skipping sample '{samp}' because sample IDs are not all defined.")
        continue

    # Prepare a figure with 2 rows × 2 columns
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8), sharex=True)
    # Axes arrangement:
    #   (0,0) => top-left  [2D Ve]
    #   (0,1) => top-right [3D Ve]
    #   (1,0) => bot-left  [2D 2nd diff]
    #   (1,1) => bot-right [3D 2nd diff]

    # For coloring by z-step, let's pick a colormap or palette
    # (or use line styles, etc.)
    colors = sns.color_palette("viridis", n_colors=len(z_steps))

    for cidx, zstep in enumerate(z_steps):
        df_ve, df_2diff = get_alldata_for_zsteps(
            sample_name=samp,
            samp2d_ud=ids['2d_ud'],
            samp3d_ud=ids['3d_ud'],
            samp2d_def=ids['2d_def'],
            samp3d_def=ids['3d_def'],
            zstep=zstep,
        )
        # df_ve has columns: [z, fiber, ve, dim, deformed]
        # df_2diff has columns: [z, fiber, ve, dim, deformed]

        # -------------------- PLOT Ve -------------------------
        # 2D => top-left
        ax_2d_ve = axes[0, 0]
        sns.lineplot(
            data=df_ve[df_ve['dim'] == '2d'],
            x='z',
            y='ve',
            style='deformed',
            errorbar=('sd', 1),
            ax=ax_2d_ve,
            color=colors[cidx],
            label=f'zstep={zstep}' if cidx == 0 else "",
        )
        ax_2d_ve.set_title(f"{samp}: 2D Ve")

        # 3D => top-right
        ax_3d_ve = axes[0, 1]
        sns.lineplot(
            data=df_ve[df_ve['dim'] == '3d'],
            x='z',
            y='ve',
            style='deformed',
            errorbar=('sd', 1),
            ax=ax_3d_ve,
            color=colors[cidx],
            label=f'zstep={zstep}' if cidx == 0 else "",
        )
        ax_3d_ve.set_title(f"{samp}: 3D Ve")

        # -------------------- PLOT 2diff -------------------------
        # 2D => bottom-left
        ax_2d_2diff = axes[1, 0]
        sns.lineplot(
            data=df_2diff[df_2diff['dim'] == '2d'],
            x='z',
            y='ve',
            style='deformed',
            errorbar=('sd', 1),
            ax=ax_2d_2diff,
            color=colors[cidx],
            label=f'zstep={zstep}' if cidx == 0 else "",
        )
        ax_2d_2diff.set_title(f"{samp}: 2D 2nd diff")

        # 3D => bottom-right
        ax_3d_2diff = axes[1, 1]
        sns.lineplot(
            data=df_2diff[df_2diff['dim'] == '3d'],
            x='z',
            y='ve',
            style='deformed',
            errorbar=('sd', 1),
            ax=ax_3d_2diff,
            color=colors[cidx],
            label=f'zstep={zstep}' if cidx == 0 else "",
        )
        ax_3d_2diff.set_title(f"{samp}: 3D 2nd diff")

    # Label axes
    axes[1, 0].set_xlabel('z (μm)')
    axes[1, 1].set_xlabel('z (μm)')
    axes[0, 0].set_ylabel('Ve (mV)')
    axes[1, 0].set_ylabel(r'$d^2Ve/dx^2$ (mV/μm$^2$)')

    # Legend handling (put them in one place, or keep separate)
    # The `style='deformed'` from Seaborn typically tries to create separate legend entries,
    # plus we have color by z-step.  One approach is to unify them:
    for axrow in axes:
        for ax in axrow:
            ax.legend().set_visible(False)

    # Combine color + style legends in a single location
    # You can build a custom legend if you want
    fig.legend(loc='upper right', bbox_to_anchor=(1.1, 0.95))

    plt.tight_layout()
    plt.show()
