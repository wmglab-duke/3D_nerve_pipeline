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

# diams = [3, 13]

# n_sims = [0, 5]

# TODO double check that the ordering of everything is correct


samp3d = 673
samplename = '6R'
n_sim = 0
simnum = 3
threshdata = pd.concat(
    [pd.read_csv(f'thresh_unmatched_sim{simnum}_og.csv'), pd.read_csv(f'thresh_unmatched_sim{simnum}_def.csv')]
).query(f'nerve_label=="{samplename}"')
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
dz = 300 if n_sim == 0 else 1300

for mfi, value in data.items():
    ves = []
    contacts = [0, 1]
    weights = [1, -1] if simnum == 3 else [0, -1]
    for contact, weight in zip(contacts, weights):
        threedpath = rf'D:\threed_ascent\samples\{samp3d}\models\0\sims\3\ss_bases\{contact}\{mfi}.dat'

        ves.append(np.loadtxt(threedpath, skiprows=1) * 1000 * weight)  # convert to mv and cathodic
    value['ve'] = np.sum(ves, axis=0)  # TODO should use the actual input ve to account for fiber jitter

    threedc = rf'D:\threed_ascent\samples\{samp3d}\models\0\sims\3\ss_coords\{mfi}.dat'

    value['arc_coords'] = np.loadtxt(threedc, skiprows=1)[:, 2]
    # TODO use actual fiber node coords
    coordfib = rf'D:\threed_ascent\samples\{samp3d}\models\0\sims\3\fibersets\{n_sim}\{mfi}.dat'
    value['node_coords'] = np.loadtxt(coordfib, skiprows=1)[:, 2][::11]

    threedcoord = rf'D:\threed_ascent\samples\{samp3d}\models\0\sims\3\3D_fiberset\{mfi}.dat'
    value['threed_coords'] = np.loadtxt(threedcoord, skiprows=1)[:, 2]

    # new array node_ve is ve at each node coord (find nearest arc coord then get ve)
    value['node_ve'] = np.zeros(len(value['node_coords']))
    for i, node in enumerate(value['node_coords']):
        idx = (np.abs(value['arc_coords'] - node)).argmin()
        value['node_ve'][i] = value['ve'][idx]

    # now get z coords for each node
    value['node_z_coords'] = np.zeros(len(value['node_coords']))
    for i, node in enumerate(value['node_coords']):
        idx = (np.abs(value['arc_coords'] - node)).argmin()
        value['node_z_coords'][i] = value['z_coords'][idx]

    # downsample to every nth point
    value['ve'] = value['ve'][::dz]
    value['arc_coords'] = value['arc_coords'][::dz]
    value['z_coords'] = value['z_coords'][::dz]
    # value['3d_coords'] = value['3d_coords'][::dz]
    # calculate second derivative
    value['ve_2diff'] = np.diff(np.diff(value['ve']))
    value['2diff_z_coords'] = value['z_coords'][1:-1]
    # now for node ve
    value['node_ve_2diff'] = np.diff(np.diff(value['node_ve']))
    value['node_2diff_z_coords'] = value['node_z_coords'][1:-1]
    # find z position of max node_ve_2diff and min node_ve_2diff
    value['max_node_2diff_z'] = value['node_2diff_z_coords'][np.argmax(value['node_ve_2diff'])]
    value['min_node_2diff_z'] = value['node_2diff_z_coords'][np.argmin(value['node_ve_2diff'])]

# %% get the first difference of peri thk
print('slides')
# loop through each slide, and calculate z position
for i, slide in enumerate(slidelist):
    zpos = i * 20  # 20 microns per slice
    if zpos < dz * 2 or zpos > 50100 - dz * 2:
        continue  # use this to check only a specific range
    if zpos < 15000 or zpos > 35000:
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
# specmfi = [random.randint(0,200) for i in range(10)]
# specmfi=[]
specmfi = [176, 177]
for name, cuffspan in zip(['approach', 'recede'], [[25000, 29000], [29000, 34000]]):
    lr = int(2000 / dz)
    threshthk = 25000
    summed_2diffs_split = np.zeros(lr * 2)
    summed_2diffs_merge = np.zeros(lr * 2)
    normalizing = True
    absing = False
    # fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    mfis = [random.randint(0, 200) for i in range(10)]
    actfig, actax = plt.subplots()
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
            idx = (np.abs(value["node_2diff_z_coords"] - merge)).argmin()
            vals = value["node_ve_2diff"][idx - lr : idx + lr]
            if normalizing:
                vals = vals / np.amax(np.abs(vals))
            if absing:
                vals = np.abs(vals)
            summed_2diffs_merge += vals
            # axs[0].plot(np.arange(-lr * dz, lr * dz, dz), vals, 'k')
        for split in split_locs:
            idx = (np.abs(value["node_2diff_z_coords"] - split)).argmin()
            vals = value["node_ve_2diff"][idx - lr : idx + lr]
            if normalizing:
                vals = vals / np.amax(np.abs(vals))
            if absing:
                vals = np.abs(vals)
            summed_2diffs_split += vals
            # axs[1].plot(np.arange(-lr * dz, lr * dz, dz), vals, 'k')
        if mfi in specmfi:
            figfib, axfib = plt.subplots(2, 1, sharex=True, sharey=False)
            axfib[0].plot(value["node_2diff_z_coords"], value["node_ve_2diff"], 'k', label='3D', marker='o')
            # plt.plot(inner_coords + ds / 2, diff2ds, 'r--', label='2D')
            axfib[0].axvspan(28000, 30000, color='r', alpha=0.2, label='contact')
            axfib[0].axvspan(26800, 31200, color='k', alpha=0.2, label='insulation')
            axfib[0].axvspan(20000, 22000, color='r', alpha=0.2, label='_contact')
            axfib[0].axvspan(18800, 23200, color='k', alpha=0.2, label='_insulation')

            figfib.set_size_inches(12, 8)
            axfib[1].plot(value["inner_coords"], value["inner_areas"], 'b', label='area', alpha=0.65)
            axfib[1].set_ylabel('fascicle area', color='b')
            axfib[1].tick_params(axis='y', labelcolor='b')
            axfib[0].set_ylabel(f'Second difference of Ve\nfiber diam = {3 if n_sim==0 else 13}', color='k')
            axfib[0].set_title(f'Fiber {mfi}, sample {samp3d}')
            axfib[0].axhline(0, color='gray', linestyle=':')
            # add vertical lines for merges and splits
            for merge in merge_locs_raw:
                axfib[0].axvline(merge, color='r', linestyle='--', label='merge' if merge == merge_locs_raw[0] else '_')
            for split in split_locs_raw:
                axfib[0].axvline(split, color='k', linestyle='--', label='split' if split == split_locs_raw[0] else '_')
            plt.xlim(15000, 35000)
            axfib[0].legend(loc='lower right')
            # obtain and plot ap init location
            for nsim in [n_sim]:
                init_loc = threshdata.query(f'master_fiber_index == {mfi} and sample == {samp3d} and nsim=={nsim}')
                assert len(init_loc) == 1
                zpos = init_loc['activation_zpos'].values[0]
                # add star to second diff plot, get y value from finding nearest value in z coords, then take y value from ve_2diff
                yval = value["node_ve_2diff"][(np.abs(value["node_2diff_z_coords"] - zpos)).argmin()]
                axfib[0].plot(zpos, yval, 'r*', markersize=10, label='activation location')
        for nsim in [n_sim]:
            init_loc = threshdata.query(f'master_fiber_index == {mfi} and sample == {samp3d} and nsim=={nsim}')
            assert len(init_loc) == 1
            zpos = init_loc['activation_zpos'].values[0]
            # add star to second diff plot, get y value from finding nearest value in z coords, then take y value from ve_2diff
            yval = value["node_ve_2diff"][(np.abs(value["node_2diff_z_coords"] - zpos)).argmin()]
            actax.plot(zpos, yval, 'k*', markersize=10)
    actax.set_ylabel("Second Difference")
    actax.set_xlabel("Initiation z-position (μm)")
    actax.axvspan(28000, 30000, color='r', alpha=0.2, label='contact')
    actax.axvspan(26800, 31200, color='k', alpha=0.2, label='insulation')
    actax.axvspan(20000, 22000, color='r', alpha=0.2, label='_contact')
    actax.axvspan(18800, 23200, color='k', alpha=0.2, label='_insulation')
    actax.axhline(0, color='gray', linestyle=':')
    actax.legend()
    actfig.show()
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
plt.axvspan(25000, 29000, color='b', alpha=0.2, label='approach')
plt.axvspan(29000, 34000, color='orange', alpha=0.2, label='recede')
plt.legend()

# %%
import seaborn as sns

# sns.catplot(data=threshdata, kind='strip', col='nsim', y='activation_zpos')#%% now compile and plot all peak and min 2nd diffs as histograms
# first, compile all peak and min 2nd diffs
allpeaks = [v['max_node_2diff_z'] for v in data.values()]
allmins = [v['min_node_2diff_z'] for v in data.values()]
# plot using seaborn histogram
sns.histplot(
    allpeaks,
    fill=True,
    alpha=0.1,
    label='+peaks',
    element='step',
)
sns.histplot(
    allmins,
    fill=True,
    alpha=0.1,
    label='-peaks',
    element='step',
)
plt.axvspan(28000, 30000, color='r', alpha=0.2, label='contact')
plt.axvspan(26800, 31200, color='k', alpha=0.2, label='insulation')
plt.axvspan(20000, 22000, color='r', alpha=0.2, label='_contact')
plt.axvspan(18800, 23200, color='k', alpha=0.2, label='_insulation')
plt.legend()
plt.xlim(18000, 32000)
