import os
import pickle
import sys

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from shapely.geometry import Point

font = {'size': 12}

matplotlib.rc('font', **font)

sys.path.append(r'C:\nrn\lib\python')  # noqa: E800

os.chdir('../..')

import pandas as pd
from shapely.strtree import STRtree

threshdata = pd.concat([pd.read_csv('thresh_unmatched_sim3_og.csv'), pd.read_csv('thresh_unmatched_sim3_def.csv')])

# diams = [3, 13]
allsum, alldiff = [], []
# n_sims = [0, 5]
simnum = 3
threeds = [2531, 3731, 5731, 6731]
samps = ['2Ldef', '3Rdef', '5Rdef', '6Rdef']
# threeds = [253,273,373,573,653,673]
# samps = ['2L', '2R', '3R', '5R','6L','6R']
names, spans = ['approach1', 'recede1', 'approach2', 'recede2'], [
    [10000, 21050],
    [21050, 25000],
    [25000, 29050],
    [29050, 40000],
]
# specmfi = [random.randint(0, 200) for i in range(10)]
specmfi = [26]
normalizing = True
absing = False
for samp3d, samplename in zip(threeds, samps):
    # %% load fiber data
    print('dataload')
    fiberpath = os.path.join(os.getcwd(), fr'samples\{samp3d}\models\0\sims\3\3D_fiberset')
    slidespath = os.path.join(os.getcwd(), fr'input\slides\{samplename}slides.obj')

    # load pickled slidelist
    with open(slidespath, 'rb') as f:
        slidelist = pickle.load(f)

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
    dz = 20

    for mfi, value in data.items():
        ves = []
        contacts = [0, 1]
        weights = [1, -1] if simnum == 3 else [0, -1]
        for contact, weight in zip(contacts, weights):
            threedpath = rf'D:\threed_final\samples\{samp3d}\models\0\sims\3\ss_bases\{contact}\{mfi}.dat'

            ves.append(np.loadtxt(threedpath, skiprows=1) * 1000 * weight)  # convert to mv and cathodic
        value['ve'] = np.sum(ves, axis=0)  # TODO should use the actual input ve to account for fiber jitter

        threedc = rf'D:\threed_final\samples\{samp3d}\models\0\sims\3\ss_coords\{mfi}.dat'

        value['arc_coords'] = np.loadtxt(threedc, skiprows=1)[:, 2]

        # downsample to every nth point
        value['ve'] = value['ve'][::dz]
        value['arc_coords'] = value['arc_coords'][::dz]
        value['z_coords'] = value['z_coords'][::dz]
        # value['3d_coords'] = value['3d_coords'][::dz]
        # calculate second derivative
        value['ve_2diff'] = np.diff(np.diff(value['ve']))
        value['2diff_z_coords'] = value['z_coords'][1:-1]
    # %% get the first difference of peri thk
    print('slides')
    # loop through each slide, and calculate z position
    for i, slide in enumerate(slidelist):
        zpos = i * 20  # 20 microns per slice
        if zpos < dz * 2 or zpos > 50100 - dz * 2:
            continue  # use this to check only a specific range
        # if zpos < 20000 or zpos > 35000:
        #     continue  # use this to check only a specific range
        # get list of x,y coordinates for each fiber at this z position
        xs = []
        ys = []
        for mfi, value in data.items():
            fiber = value['3d_coords']
            # find the index of the closest z value
            idx = (np.abs(fiber[:, 2] - zpos)).argmin()
            # get the x,y coordinates at that index
            x = fiber[idx, 0]
            y = fiber[idx, 1]
            # append to slide list
            xs.append(x)
            ys.append(-y)
            # find which inner fascicle is in
            for inner in slide.inners():
                if inner.contains(Point(x, -y)):
                    value["inner_areas"].append(inner.area())
                    # inner_2diffs.append(ve3_2diff[idx // ds]) #TODO
                    value["inner_coords"].append(zpos)
                    break
            else:
                # try finding nearest inner using STRtree
                tree = STRtree(slide.inners(polygon=True))
                nearest = tree.nearest(Point(x, -y))
                if nearest.distance(Point(x, -y)) < 50:
                    value["inner_areas"].append(nearest.area)
                    # inner_2diffs.append(ve3_2diff[idx // ds])
                    value["inner_coords"].append(zpos)
                else:
                    plt.plot(x, -y, 'rx')
                    slide.plot()
                    raise ValueError('No inner found')
            value["inner_coords"] = value["inner_coords"]
    # %% now do cross correlelograms
    print('plotgen')
    for name, cuffspan in zip(names, spans):
        lr = int(1000 / dz)
        threshthk = 25000
        summed_2diffs_split = np.zeros(lr * 2)
        summed_2diffs_merge = np.zeros(lr * 2)
        for mfi, value in data.items():
            # calculate first difference of inner areas
            value["inner_1diffs"] = np.diff(value["inner_areas"])
            # calculate new inner coords
            value["inner_1diff_coords"] = np.array(value["inner_coords"])[1:] - 10  # half of z step
            # first, threshold the absolute value of perineurium thickness to find splits and merges
            merges = np.where(value["inner_1diffs"] > threshthk)[0]
            splits = np.where(value["inner_1diffs"] < -threshthk)[0]
            # replae with find peaks
            merges = find_peaks(value["inner_1diffs"], height=threshthk, distance=5)[0]
            splits = find_peaks(-value["inner_1diffs"], height=threshthk, distance=5)[0]
            # get the x locations of the merges and splits
            merge_locs_raw = value["inner_1diff_coords"][merges]
            split_locs_raw = value["inner_1diff_coords"][splits]
            # only consider merges and spits from 20000 to 40000
            merge_locs = merge_locs_raw[(merge_locs_raw > cuffspan[0]) & (merge_locs_raw < cuffspan[1])]
            split_locs = split_locs_raw[(split_locs_raw > cuffspan[0]) & (split_locs_raw < cuffspan[1])]

            # find the closest fiber point to each merge, and get the second difference to the left and right using lr
            for merge in merge_locs:
                idx = (np.abs(value["2diff_z_coords"] - merge)).argmin()
                vals = value["ve_2diff"][idx - lr : idx + lr]
                if normalizing:
                    vals = vals / np.amax(np.abs(vals))
                if absing:
                    vals = np.abs(vals)
                summed_2diffs_merge += vals
            for split in split_locs:
                idx = (np.abs(value["2diff_z_coords"] - split)).argmin()
                vals = value["ve_2diff"][idx - lr : idx + lr]
                if normalizing:
                    vals = vals / np.amax(np.abs(vals))
                if absing:
                    vals = np.abs(vals)
                summed_2diffs_split += vals
            if mfi in specmfi and name == names[0]:
                figfib, axfib = plt.subplots(2, 1, sharex=True, sharey=False)
                axfib[0].plot(value["2diff_z_coords"], value["ve_2diff"], 'k', label='3D')
                # plt.plot(inner_coords + ds / 2, diff2ds, 'r--', label='2D')
                axfib[0].axvspan(28000, 30000, color='r', alpha=0.2, label='contact')
                axfib[0].axvspan(26800, 31200, color='k', alpha=0.2, label='insulation')
                axfib[0].axvspan(20000, 22000, color='r', alpha=0.2, label='_contact')
                axfib[0].axvspan(18800, 23200, color='k', alpha=0.2, label='_insulation')
                figfib.set_size_inches(12, 8)
                axfib[1].plot(value["inner_coords"], value["inner_areas"], 'b', label='area', alpha=0.65)
                axfib[1].set_ylabel('fascicle area (um^2)')
                axfib[1].set_xlabel('Distance along nerve (μm)')
                axfib[0].set_ylabel(f'Second difference of Ve\ndz = {dz} μm', color='k')
                axfib[0].set_title(f'Fiber {mfi} - sample {samp3d}')
                # add vertical lines for merges and splits
                for merge in merge_locs_raw:
                    axfib[0].axvline(
                        merge, color='r', linestyle='--', label='merge' if merge == merge_locs_raw[0] else '_'
                    )
                for split in split_locs_raw:
                    axfib[0].axvline(
                        split, color='k', linestyle='--', label='split' if split == split_locs_raw[0] else '_'
                    )
                # plt.xlim(15000, 35000)
                axfib[0].legend(loc='lower right')
                axfib[0].axhline(0, color='gray', linestyle=':')
        allsum.append(summed_2diffs_split)
        alldiff.append(summed_2diffs_merge)
# %% plot each summed 2diff
fig, ax = plt.subplots(1, len(names), sharex=True, sharey=False)
# each summed 2diff is in order of nerve, then approach/recede (so 0 is approach for one nerve, 1 is recede for that nerve, 2 is approach for next nerve, etc)
# plot all approaches and recedes for sum and diff on two plots, one for approach and one for recede
for i in range(0, len(allsum), len(names)):
    for j in range(len(names)):
        ax[j].plot(np.arange(-lr * dz, lr * dz, dz), allsum[i + j], 'k', label=names[j] if i == 0 else '_')
        ax[j].plot(np.arange(-lr * dz, lr * dz, dz), alldiff[i + j], 'r', label=names[j] if i == 0 else '_')
        plt.xlabel('z position (microns)')
        plt.ylabel('Second difference of Ve')
        ax[j].axvline(x=0, color='k', linestyle='--')
        ax[j].set_title(f'{names[j]}, dz={dz}')
        plt.legend()
fig.set_size_inches(24, 4)
# %% specific plot
from itertools import cycle

cycol = cycle('bgrcmk')
plt.figure()
# plot vm of fiber 0
plt.plot(data[0]["z_coords"], data[0]["ve"] / np.abs(np.amax(data[0]["ve"])), 'k', label='3D Ve')
plt.xlabel('z position (microns)')
plt.ylabel('Normalized Value')
plt.title('Vm of fiber 0')

# plt.gca().twinx().plot(data[0]["z_coords"][1:-1], np.diff(data[0]["ve"],n=2), 'gray')
twoodvpath = rf'D:\threed_final\samples\672\models\0\sims\{simnum}\n_sims\0\data\inputs\inner0_fiber0.dat'
twoodv = -np.loadtxt(twoodvpath)[1:][::11]
twoodv_1 = np.diff(twoodv)
twoodv_2 = np.diff(twoodv, n=2)

zcoord = np.loadtxt(r'D:\threed_final\samples\672\models\0\sims\3\fibersets\0\0.dat', skiprows=1)[:, 2][::11]
plt.plot(zcoord, twoodv / np.amax(np.abs(twoodv)), label='2D Ve')
plt.plot(zcoord[1:], twoodv_1 / np.amax(np.abs(twoodv_1)), label='2D Ve 1st diff')
plt.plot(zcoord[1:-1], twoodv_2 / np.amax(np.abs(twoodv_2)), label='2D Ve 2nd diff')

# shade cuffspans
for name, span in zip(names, spans):
    plt.axvspan(span[0], span[1], color=next(cycol), alpha=0.2, label=name)
plt.legend(bbox_to_anchor=(1, 1))
