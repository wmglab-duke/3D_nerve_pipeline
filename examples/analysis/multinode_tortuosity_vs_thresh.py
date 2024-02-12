import os
import pickle
import random
import sys

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from shapely.geometry import Point

sys.path.append(r'C:\nrn\lib\python')  # noqa: E800

os.chdir('../..')

import pandas as pd
from shapely.strtree import STRtree

threshdata = pd.concat([pd.read_csv('thresh_unmatched_sim3_og.csv')]).query('type=="3D"')
tort_data = []
r2s = []
for samp3d, samplename in zip([253, 273, 373, 573, 653, 673], ['2L', '2R', '3R', '5R', '6L', '6R']):
    simnum = 3

    # %% load fiber data
    print('dataload')
    fiberpath = os.path.join(os.getcwd(), fr'samples\{samp3d}\models\0\sims\3\3D_fiberset')

    fibers = {}
    iters = 0
    # load each fiber file and append to list
    for file in os.listdir(fiberpath):
        # if file!='212.dat':
        #     continue
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

    for mfi, value in data.items():
        print(mfi)
        ves = []

        threedc = rf'D:\threed_ascent\samples\{samp3d}\models\0\sims\3\ss_coords\{mfi}.dat'

        value['arc_coords'] = np.loadtxt(threedc, skiprows=1)[:, 2]

        for n_sim in [0, 5]:
            coordfib = rf'D:\threed_ascent\samples\{samp3d}\models\0\sims\3\fibersets\{n_sim}\{mfi}.dat'
            value[f'node_coords_{n_sim}'] = np.loadtxt(coordfib, skiprows=1)[:, 2][::11]

            threedcoord = rf'D:\threed_ascent\samples\{samp3d}\models\0\sims\3\3D_fiberset\{mfi}.dat'
            value[f'threed_coords_{n_sim}'] = np.loadtxt(threedcoord, skiprows=1)[:, 2]

            value[f'node_3d_coords_{n_sim}'] = np.zeros([len(value[f'node_coords_{n_sim}']), 3])
            for i, node in enumerate(value[f'node_coords_{n_sim}']):
                idx = (np.abs(value['arc_coords'] - node)).argmin()
                value[f'node_3d_coords_{n_sim}'][i] = value['3d_coords'][idx]

    # %% now calculate various tortuosities
    print('tortuosities')
    for mfi, value in data.items():
        for n_sim in [0, 5]:
            # calculate the tortuosity of the nodes around where the action potential initated
            # first get the data for this fiber
            fiberdat = threshdata.query(f'master_fiber_index=={mfi} and nsim=={n_sim} and sample=={samp3d}')
            assert len(fiberdat) == 1
            apnode = fiberdat['apnode'].values[0]
            thresh = fiberdat['threshold'].values[0]
            for nodeno in [1, 2, 3, 4, 5]:
                # get the node coordinates around the AP node
                nodecoords = value[f'node_3d_coords_{n_sim}'][apnode - nodeno : apnode + nodeno + 1]
                # tortuosity is arc length /  end to end distance
                tortuosity = np.sum(np.linalg.norm(np.diff(nodecoords, axis=0), axis=1)) / np.linalg.norm(
                    nodecoords[-1] - nodecoords[0]
                )
                tort_data.append(
                    {
                        'tortuosity': tortuosity,
                        'nodeno': nodeno,
                        'mfi': mfi,
                        'n_sim': n_sim,
                        'thresh': thresh,
                        'sample': samplename,
                    }
                )
            # now do for the tortuosity of the whole fiber
            tortuosity = np.sum(
                np.linalg.norm(np.diff(value[f'node_3d_coords_{n_sim}'], axis=0), axis=1)
            ) / np.linalg.norm(value[f'node_3d_coords_{n_sim}'][-1] - value[f'node_3d_coords_{n_sim}'][0])
            tort_data.append(
                {
                    'tortuosity': tortuosity,
                    'nodeno': float('Inf'),
                    'mfi': mfi,
                    'n_sim': n_sim,
                    'thresh': thresh,
                    'sample': samplename,
                }
            )
# %% plot
tort_data = pd.DataFrame(tort_data)
# plot tortuosity vs threshold for each nodeno
import seaborn as sns

sns.set(font_scale=1.75, style='whitegrid')
g = sns.lmplot(
    x='tortuosity',
    y='thresh',
    data=tort_data,
    row='n_sim',
    sharex=False,
    sharey=False,
    col='nodeno',
    facet_kws=dict(margin_titles=True),
    hue='sample',
    palette='colorblind',
)
g.set_titles(col_template='Node distance {col_name}', row_template='')
g.axes[0, 0].set_ylabel('Threshold (mA)\nD=3 μm')
g.axes[1, 0].set_ylabel('Threshold (mA)\nD=13 μm')
# %% calculate and plot R2 values between each tortuosity and threshold
from scipy.stats import linregress

r2s = []
for nodeno in [1, 2, 3, 4, 5, float('Inf')]:
    for n_sim in [0, 5]:
        for samplename in ['2L', '2R', '3R', '5R', '6L', '6R']:
            torts = tort_data.query(f'n_sim=={n_sim} and nodeno=={nodeno} and sample=="{samplename}"')['tortuosity']
            threshs = tort_data.query(f'n_sim=={n_sim} and nodeno=={nodeno} and sample=="{samplename}"')['thresh']
            r2s.append(
                {'r2': linregress(torts, threshs)[2] ** 2, 'nodeno': nodeno, 'n_sim': n_sim, 'sample': samplename}
            )
r2s = pd.DataFrame(r2s)
sns.set(font_scale=1.75, style='whitegrid')
# use lineplot
for n_sim in [0, 5]:
    plt.figure()
    g = sns.pointplot(x='nodeno', y='r2', data=r2s.query(f'n_sim=={n_sim}'), palette='colorblind', hue='sample')
    g.set_xlabel('Node distance from AP node')
    g.set_ylabel('R$^2$')
    sns.move_legend(g, [1, 0])
    plt.title(n_sim)
    # change legend text to 3 and 13
    # handles, labels = g.get_legend_handles_labels()
    # g.legend(handles, ['3','13'],title='D (μm)')
