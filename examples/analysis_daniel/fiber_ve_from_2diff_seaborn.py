import os
import sys

import matplotlib
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

font = {'size': 12}

matplotlib.rc('font', **font)

sys.path.append(r'C:\nrn\lib\python')  # noqa: E800

os.chdir('../..')

import pandas as pd

threshdata = pd.concat([pd.read_csv('thresh_unmatched_sim10_og.csv'), pd.read_csv('thresh_unmatched_sim10_def.csv')])
pal2d3d = ['#d95f02', '#7570b3']

samp2d_ud = 571
samp3d_ud = 573
samplename = '5R'
samp2ddef = 5711
samp3ddef = 5731
alldata = {}
for samp2d, samp3d, deformed in zip([samp2d_ud, samp2ddef], [samp3d_ud, samp3ddef], ['og', 'def']):
    # %% load fiber data
    print('dataload')
    fiberpath = os.path.join(os.getcwd(), fr'samples\{samp3d}\models\0\sims\3\3D_fiberset')
    slidespath = os.path.join(os.getcwd(), fr'input\slides\{samplename}slides.obj')

    fibers = {}
    iters = 0
    # load each fiber file and append to list
    for file in os.listdir(fiberpath):
        if file.endswith('.dat'):
            fibers[int(file.replace('.dat', ''))] = np.loadtxt(os.path.join(fiberpath, file), skiprows=1)
        iters += 1
        # if iters>=50:
        #     break
    # %%
    print('preproc')
    data = {
        key: {"z_coords": value[:, 2], "3d_coords": value, "inner_areas": [], "inner_coords": []}
        for key, value in fibers.items()
    }

    # define z coords of interest
    zs = np.arange(50, 50100, 1000)

    for mfi, value in data.items():
        ves = []
        contacts = [0, 1]
        weights = [1, -1]
        for contact, weight in zip(contacts, weights):
            threedpath = rf'D:\threed_final\samples\{samp3d}\models\0\sims\3\ss_bases\{contact}\{mfi}.dat'

            ves.append(np.loadtxt(threedpath, skiprows=1) * 1000 * weight)  # convert to mv and cathodic
        value['ve'] = np.sum(ves, axis=0)  # TODO should use the actual input ve to account for fiber jitter

        threedc = rf'D:\threed_final\samples\{samp3d}\models\0\sims\3\ss_coords\{mfi}.dat'

        value['arc_coords'] = np.loadtxt(threedc, skiprows=1)[:, 2]

        threedpath = rf'D:\threed_final\samples\{samp3d}\models\0\sims\3\ss_bases\1\{mfi}.dat'

        threedc = rf'D:\threed_final\samples\{samp3d}\models\0\sims\3\ss_coords\{mfi}.dat'

        # calculate second derivative
        value['ve_2diff'] = np.diff(np.diff(value['ve']))
        value['2diff_z_coords'] = value['z_coords'][1:-1]

        # now 2d
        ves = []
        for contact, weight in zip(contacts, weights):
            twoodpath = rf'D:\threed_final\samples\{samp2d}\models\0\sims\3\ss_bases\{contact}\{mfi}.dat'

            ves.append(np.loadtxt(twoodpath, skiprows=1) * 1000 * weight)  # convert to mv and cathodic
        value['ve_2d'] = np.sum(ves, axis=0)  # TODO should use the actual input ve to account for fiber jitter

        twoodc = rf'D:\threed_final\samples\{samp2d}\models\0\sims\3\ss_coords\{mfi}.dat'

        value['z_coords_2d'] = np.loadtxt(twoodc, skiprows=1)[:, 2]

        # calculate second derivative
        value['ve_2diff_2d'] = np.diff(np.diff(value['ve_2d']))
        value['2diff_z_coords_2d'] = value['z_coords_2d'][1:-1]

        # for both 2d and 3d, find the indices of the z_coords at the zs of interest, and get the ve at those indices
        ve_2_inds = [np.argmin(np.abs(value['z_coords_2d'] - z)) for z in zs]
        value['ve_2d_at_zs'] = value['ve_2d'][ve_2_inds]
        ve_3_inds = [np.argmin(np.abs(value['z_coords'] - z)) for z in zs]
        value['ve_3d_at_zs'] = value['ve'][ve_3_inds]
    # create two wide form dfs, one for 2d and one for 3d. Each row is a fiber, each column is a z coord, and the value is the ve at that z coord
    df_2d = pd.DataFrame({key: value['ve_2d_at_zs'] for key, value in data.items()})
    df_3d = pd.DataFrame({key: value['ve_3d_at_zs'] for key, value in data.items()})
    # set the column names to the z coords
    df_2d.index = zs
    df_3d.index = zs
    # convert to long form
    df_2d = df_2d.reset_index().melt(id_vars='index')
    df_3d = df_3d.reset_index().melt(id_vars='index')
    # rename columns
    df_2d.columns = ['z', 'fiber', 've']
    df_3d.columns = ['z', 'fiber', 've']
    # add a column for 2d or 3d
    df_2d['dim'] = '2d'
    df_3d['dim'] = '3d'
    # concatenate
    df = pd.concat([df_2d, df_3d])

    alldata[deformed] = df
# %%
print('plot')
# convert alldata to long form dataframe and add a column for deformed or not
df = pd.concat([value.assign(deformed=key) for key, value in alldata.items()])
df['deformed'] = df.deformed.replace({"og": 'no', "def": 'yes'})
# plot seaborn
g = sns.lineplot(
    data=df, x='z', y='ve', hue='dim', style='deformed', err_style='band', errorbar=('sd', 1), palette=pal2d3d
)
sns.move_legend(g, loc='upper left')
plt.xlabel('z (μm)')
plt.ylabel('Ve (mV)')
plt.title(samplename)
sns.relplot(data=df, x='z', y='ve', hue='dim', col='deformed', kind='line', estimator=None, units='fiber')
