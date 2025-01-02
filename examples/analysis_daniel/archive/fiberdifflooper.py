import sys

import numpy as np
from shapely.geometry import Point

sys.path.append(r'C:\nrn\lib\python')  # noqa: E800
import os

os.chdir('../..')


sim = 3

# diams = [3, 13]

# n_sims = [0, 5]


sample = 57

innernum, fibernum, mfi = 11, 1, 26

mfis = [0, 1, 2, 3, 4]

diams = [3] * len(mfis)

n_sims = [0] * len(mfis)

samplename = '5R'
for dz in [1]:
    for d, n, mfi in zip(diams, n_sims, mfis):
        # twood = (
        #     rf'D:\threed_ascent\samples\{sample}2\models\0\sims\{sim}\n_sims\{n}\data\inputs\inner{innernum}_fiber{fibernum}.dat'
        # )
        threed = rf'D:\threed_ascent\samples\{sample}3\models\0\sims\{sim}\n_sims\{n}\data\inputs\inner0_fiber{mfi}.dat'

        # plot both fibers ve

        # ve2 = np.loadtxt(twood, skiprows=1)[::11] * 1000
        ve3 = np.loadtxt(threed, skiprows=1)[::-1][::11] * 1000

        import matplotlib.pyplot as plt

        # now plot the basis
        # twood = rf'D:\threed_ascent\samples\{sample}2\models\0\sims\3\ss_bases\1\{mfi}.dat'
        threed = rf'D:\threed_ascent\samples\{sample}3\models\0\sims\3\ss_bases\1\{mfi}.dat'

        # ve2 = np.loadtxt(twood, skiprows=1) * 1000
        ve3 = np.loadtxt(threed, skiprows=1)[::-1] * 1000

        twoodc = rf'D:\threed_ascent\samples\{sample}2\models\0\sims\3\ss_coords\{mfi}.dat'
        threedc = rf'D:\threed_ascent\samples\{sample}3\models\0\sims\3\ss_coords\{mfi}.dat'

        coords2 = np.loadtxt(twoodc, skiprows=1)
        coords3 = np.loadtxt(threedc, skiprows=1)

        # downsample to every nth point
        ds = dz
        # ve2 = ve2[::ds]
        ve3 = ve3[::ds]
        ve3_1diff = (np.diff(ve3))[::-1]
        ve3_2diff = (np.diff(np.diff(ve3)) * 1000)[::-1]

        # ve2_1diff = np.diff(ve2)
        # ve2_2diff = np.diff(np.diff(ve2)) * 1000
        ve3 = ve3[::-1]

        # %%
        plotdiam = d
        import os
        import pickle

        samp3d = int(str(sample) + '3')
        diamdir = rf'D:\threed_ascent\plots\testfollow\{d}'
        os.makedirs(diamdir, exist_ok=True)
        fiberpath = os.path.join(os.getcwd(), fr'samples\{samp3d}\models\0\sims\3\3D_fiberset')
        slidespath = os.path.join(os.getcwd(), fr'input\slides\{samplename}slides.obj')

        # load pickled slidelist
        with open(slidespath, 'rb') as f:
            slidelist = pickle.load(f)

        fibers = {}
        # load each fiber file and append to list
        for file in os.listdir(fiberpath):
            # print(file)
            if file != f'{mfi}.dat':
                continue
            if file.endswith('.dat'):
                fibers[int(file.replace('.dat', ''))] = np.loadtxt(os.path.join(fiberpath, file), skiprows=1)
        # %%
        inner_areas = []
        inner_2diffs = []
        inner_coords = []
        # diff2ds = []
        # loop through each slide, and calculate z position
        for i, slide in enumerate(slidelist):
            # if i % 5 != 0:  # every 100 microns
            #     continue
            zpos = i * 20  # 20 microns per slice
            # if zpos < 27000 or zpos > 31000:
            # if zpos < 25000 or zpos > 33000:
            if zpos < ds * 2 or zpos > 50100 - ds * 2:
                continue  # use this to check only a specific range
            # get list of x,y coordinates for each fiber at this z position
            xs = []
            ys = []
            for mfi, fiber in fibers.items():
                # find the index of the closest z value
                idx = (np.abs(fiber[:, 2] - zpos)).argmin()
                # get the x,y coordinates at that index
                x = fiber[idx, 0]
                y = fiber[idx, 1]
                # append to slide list
                xs.append(x)
                ys.append(-y)
                # also find index of closest z value in 2D fiber
                idx2 = (np.abs(coords2[:, -1] - zpos)).argmin()
                # find which inner fascicle is in
                for inner in slide.inners():
                    if inner.contains(Point(x, -y)):
                        inner_areas.append(inner.area())
                        inner_2diffs.append(ve3_2diff[idx // ds])
                        inner_coords.append(zpos)
                        # diff2ds.append(ve2_2diff[idx2 // ds])
                        break
                else:
                    # try finding nearest inner using STRtree
                    from shapely.strtree import STRtree

                    tree = STRtree(slide.inners(polygon=True))
                    nearest = tree.nearest(Point(x, -y))
                    if nearest.distance(Point(x, -y)) < 50:
                        inner_areas.append(nearest.area)
                        inner_2diffs.append(ve3_2diff[idx // ds])
                        inner_coords.append(zpos)
                        # diff2ds.append(ve2_2diff[idx2 // ds])
                    else:
                        raise ValueError('No inner found')
        plt.figure()
        ax = plt.gca()
        inner_coords = np.array(inner_coords)
        plt.plot(inner_coords - ds / 2, inner_2diffs, 'k', label='3D')
        # plt.plot(inner_coords + ds / 2, diff2ds, 'r--', label='2D')
        plt.axvspan(28000, 30000, color='r', alpha=0.2, label='contact')
        plt.axvspan(26800, 31200, color='k', alpha=0.2, label='insulation')

        plt.legend(loc='lower right')
        ax2 = plt.gca().twinx()
        plt.gcf().set_size_inches(12, 4)
        ax2.plot(inner_coords[1:], np.diff(inner_areas), 'b', label='area', alpha=0.65)
        ax2.set_ylabel('First difference of fascicle area', color='b')
        ax2.tick_params(axis='y', labelcolor='b')
        ax.set_ylabel('Second difference of Ve', color='k')
        plt.title(f'Fiber {mfi}')
        plt.xlim()
        # %% now do cross correlelograms
        # first, threshold the absolute value of perineurium thickness to find splits and merges

        merges = np.where(np.diff(inner_areas) > 25000)[0]
        splits = np.where(np.diff(inner_areas) < -25000)[0]
        # get the x locations of the merges and splits
        merge_locs = inner_coords[merges]
        split_locs = inner_coords[splits]
        # only consider merges and spits from 20000 to 40000
        merge_locs = merge_locs[(merge_locs > 20000) & (merge_locs < 40000)]
        split_locs = split_locs[(split_locs > 20000) & (split_locs < 40000)]
        # for each merge and split, plot ve3_2diff for 20 points before and after
        # for loc in merge_locs:
        #     idx = np.where(inner_coords == loc)[0][0]
        #     plt.figure()
        #     plt.plot(inner_coords[idx - 20:idx + 20], inner_2diffs[idx - 20:idx + 20], 'k')
        #     plt.axvline(loc, color='r')
        #     plt.title(f'Fiber {mfi} merge at {loc}')
        #     plt.xlabel('z position')
        #     plt.ylabel('Ve 2nd difference')
        # for loc in split_locs:
        #     idx = np.where(inner_coords == loc)[0][0]
        #     plt.figure()
        #     plt.plot(inner_coords[idx - 20:idx + 20], inner_2diffs[idx - 20:idx + 20], 'k')
        #     plt.axvline(loc, color='r')
        #     plt.title(f'Fiber {mfi} split at {loc}')
        #     plt.xlabel('z position')
        #     plt.ylabel('Ve 2nd difference')
        # now plot, centering on the merge/split
        plt.figure()
        for loc in merge_locs:
            idx = np.where(inner_coords == loc)[0][0]
            plt.plot(inner_coords[idx - 20 : idx + 20] - loc, inner_2diffs[idx - 20 : idx + 20], 'k')
            plt.axvline(0, color='r')
            plt.title(f'Fiber {mfi} merge at {loc}')
            plt.xlabel('z position')
            plt.ylabel('Ve 2nd difference')
        plt.figure()
        for loc in split_locs:
            idx = np.where(inner_coords == loc)[0][0]
            plt.plot(inner_coords[idx - 20 : idx + 20] - loc, inner_2diffs[idx - 20 : idx + 20], 'k')
            plt.axvline(0, color='r')
            plt.title(f'Fiber {mfi} split at {loc}')
            plt.xlabel('z position')
            plt.ylabel('Ve 2nd difference')
        plt.figure()
        # now do the same plot but sum all of the 2diffs together
        summed_2diffs = np.zeros(40)
        for loc in merge_locs:
            idx = np.where(inner_coords == loc)[0][0]
            summed_2diffs += inner_2diffs[idx - 20 : idx + 20]
        plt.plot(range(-20, 20), summed_2diffs, 'k')
        plt.figure()
        summed_2diffs = np.zeros(40)
        for loc in split_locs:
            idx = np.where(inner_coords == loc)[0][0]
            summed_2diffs += inner_2diffs[idx - 20 : idx + 20]
        plt.plot(range(-20, 20), summed_2diffs, 'r')
