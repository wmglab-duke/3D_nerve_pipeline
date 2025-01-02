import os
import pickle
import random
import sys

import matplotlib
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from shapely.geometry import Point

font = {'size': 12}

matplotlib.rc('font', **font)

sys.path.append(r'C:\nrn\lib\python')  # noqa: E800

os.chdir('../..')

import pandas as pd
from shapely.strtree import STRtree

threshdata = pd.concat([pd.read_csv('thresh_unmatched_sim333_og.csv'), pd.read_csv('thresh_unmatched_sim333_def.csv')])

# diams = [3, 13]

# n_sims = [0, 5]


samp3d = 5731
samplename = '5Rdef'
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
    threedpath = rf'D:\threed_ascent\samples\{samp3d}\models\0\sims\3\ss_bases\1\{mfi}.dat'

    value['ve'] = np.loadtxt(threedpath, skiprows=1) * -1000  # convert to mv and cathodic

    threedc = rf'D:\threed_ascent\samples\{samp3d}\models\0\sims\3\ss_coords\{mfi}.dat'

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
    if zpos < 20000 or zpos > 35000:
        continue  # use this to check only a specific range
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
for name, cuffspan in zip(['approach', 'recede'], [[9000, 29000], [29000, 49000]]):
    lr = int(2000 / dz)
    threshthk = 25000
    summed_2diffs_split = np.zeros(lr * 2)
    summed_2diffs_merge = np.zeros(lr * 2)
    normalizing = True
    absing = False
    mfis = [random.randint(0, 200) for i in range(10)]
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
        if mfi in [26]:
            figfib, axfib = plt.subplots(2, 1, sharex=True, sharey=False)
            axfib[0].plot(value["2diff_z_coords"], value["ve_2diff"], 'k', label='3D')
            # plt.plot(inner_coords + ds / 2, diff2ds, 'r--', label='2D')
            axfib[0].axvspan(28000, 30000, color='r', alpha=0.2, label='contact')
            axfib[0].axvspan(26800, 31200, color='k', alpha=0.2, label='insulation')

            figfib.set_size_inches(12, 8)
            axfib[1].plot(value["inner_coords"], value["inner_areas"], 'b', label='area', alpha=0.65)
            axfib[1].set_ylabel('fascicle area (um^2)')
            axfib[1].set_xlabel('Distance along nerve (μm)')
            axfib[0].set_ylabel(f'Second difference of Ve\ndz = {dz} μm', color='k')
            axfib[0].set_title(f'Fiber {mfi} - sample {samp3d}')
            # add vertical lines for merges and splits
            for merge in merge_locs_raw:
                axfib[0].axvline(merge, color='r', linestyle='--', label='merge' if merge == merge_locs_raw[0] else '_')
            for split in split_locs_raw:
                axfib[0].axvline(split, color='k', linestyle='--', label='split' if split == split_locs_raw[0] else '_')
            plt.xlim(20000, 35000)
            axfib[0].legend(loc='lower right')
            axfib[0].axhline(0, color='gray', linestyle=':')

    plt.figure()
    plt.plot(np.arange(-lr * dz, lr * dz, dz), summed_2diffs_split, 'k', label='split')
    plt.plot(np.arange(-lr * dz, lr * dz, dz), summed_2diffs_merge, 'r', label='merge')
    plt.xlabel('z position (microns)')
    plt.ylabel('Second difference of Ve')
    plt.axvline(x=0, color='k', linestyle='--')
    plt.title(f'dz={dz} microns, sample {samplename} - {name}')
    plt.legend()
# %% specific plot
plt.figure()
# plot vm of fiber 0
plt.plot(data[0]["z_coords"], data[0]["ve"], 'k')
plt.xlabel('z position (microns)')
plt.ylabel('Vm (mV)')
plt.title('Vm of fiber 0')
# shade cuffspans
plt.axvspan(9000, 29000, color='b', alpha=0.2, label='approach')
plt.axvspan(29000, 49000, color='orange', alpha=0.2, label='recede')
plt.legend()
