import sys

import numpy as np
from shapely.geometry import Point

sys.path.append(r'C:\nrn\lib\python')  # noqa: E800
import os

os.chdir('../..')

from moviepy.video.io.bindings import mplfig_to_npimage
from moviepy.video.VideoClip import VideoClip

sim = 3

diams = [3, 13]

n_sims = [0, 5]

sample = 67

innernum, fibernum, mfi = 0, 0, 212

samplename = '6R'

for d, n in zip(diams, n_sims):

    twood = rf'D:\threed_ascent\samples\{sample}2\models\0\sims\{sim}\n_sims\{n}\data\inputs\inner{innernum}_fiber{fibernum}.dat'
    threed = rf'D:\threed_ascent\samples\{sample}3\models\0\sims\{sim}\n_sims\{n}\data\inputs\inner0_fiber{mfi}.dat'

    # plot both fibers ve

    ve2 = -np.loadtxt(twood, skiprows=1)[::11] * 1000
    ve3 = -np.loadtxt(threed, skiprows=1)[::-1][::11] * 1000

    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(ve2, label='twood')
    plt.plot(ve3, label='threed')
    plt.legend()
    plt.title(f'Ve for {d} micron fibers')
    # note that the 3D fiber has a higher peak, at all points, not just where the perineurium thickness is smallest)

    # plot first differential
    plt.figure()
    plt.plot(np.diff(ve2), label='twood')
    plt.plot(np.diff(ve3), label='threed')
    plt.title(f'First differential of Ve for {d} micron fibers')
    plt.legend()

    # plot second differential
    plt.figure()
    plt.plot(np.diff(np.diff(ve2)), label='twood')
    plt.plot(np.diff(np.diff(ve3)), label='threed')
    plt.title(f'Second differential of Ve for {d} micron fibers')
    plt.legend()

    # now plot the basis

    twood = rf'D:\threed_ascent\samples\{sample}2\models\0\sims\3\ss_bases\1\{mfi}.dat'
    threed = rf'D:\threed_ascent\samples\{sample}3\models\0\sims\3\ss_bases\1\{mfi}.dat'

    ve2 = np.loadtxt(twood, skiprows=1) * 1000
    ve3 = np.loadtxt(threed, skiprows=1)[::-1] * 1000

    twoodc = rf'D:\threed_ascent\samples\{sample}2\models\0\sims\3\ss_coords\{mfi}.dat'
    threedc = rf'D:\threed_ascent\samples\{sample}3\models\0\sims\3\ss_coords\{mfi}.dat'

    coords2 = np.loadtxt(twoodc, skiprows=1)
    coords3 = np.loadtxt(threedc, skiprows=1)

    # downsample to every nth point
    ds = 1300
    ve2 = ve2[::ds]
    ve3 = ve3[::ds]
    ve3_1diff = (np.diff(ve3))[::-1]
    ve3_2diff = (np.diff(np.diff(ve3)) * 1000)[::-1]

    ve2_1diff = np.diff(ve2)
    ve2_2diff = np.diff(np.diff(ve2)) * 1000
    ve3 = ve3[::-1]

    #%%
    plotdiam = d
    import os
    import pickle

    samp3d = int(str(sample) + '3')
    diamdir = rf'D:\threed_ascent\plots\testfollow\{d}'
    os.makedirs(diamdir, exist_ok=True)
    fiberpath = os.path.join(os.getcwd(), r'samples\{}\models\0\sims\3\3D_fiberset'.format(samp3d))
    slidespath = os.path.join(os.getcwd(), r'input\slides\{}slides.obj'.format(samplename))

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
    diff2ds = []
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
                    diff2ds.append(ve2_2diff[idx2 // ds])
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
                    diff2ds.append(ve2_2diff[idx2 // ds])
                else:
                    raise ValueError('No inner found')
        # plot the slide and all fiber points
        if False:
            plt.figure()
            plt.scatter(xs, ys, s=6, color='k')
            # label the top of colorbar "activates first in 2D" and bottom "activates first in 3D"
            plt.title(f'Slide {i}-zpos{zpos}')
            # add colorbar with new colormap
            slide.plot(final=False)
            # plt.show()
            # state the ve at idx, round to 3 decimal places
            # print(idx,ve3[idx])
            idx2 = idx2 // ds
            idx = idx // ds
            plt.text(
                1,
                0.5,
                f'2D ve: {ve2[idx2]}\n3D ve: {ve3[idx]}',
                verticalalignment='center',
                transform=plt.gca().transAxes,
            )
            # also add text for second and first differential
            plt.text(
                1,
                0.4,
                f'2D 1st diff: {ve2_1diff[idx2]}\n3D 1st diff: {ve3_1diff[idx]}',
                verticalalignment='center',
                transform=plt.gca().transAxes,
            )
            plt.text(
                1,
                0.3,
                f'2D 2nd diff: {ve2_2diff[idx2]}\n3D 2nd diff: {ve3_2diff[idx]}',
                verticalalignment='center',
                transform=plt.gca().transAxes,
            )
            plt.savefig(os.path.join(diamdir, f'slide{i}_zpos{zpos}.png'), bbox_inches='tight')
    plt.figure()
    ax = plt.gca()
    inner_coords = np.array(inner_coords)
    plt.plot(inner_coords - ds / 2, inner_2diffs, 'k', label='3D')
    plt.plot(inner_coords + ds / 2, diff2ds, 'r--', label='2D')
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
#%%
# now simulate the 3 micron fiber 2D
from wmglab_neuron import FiberModel, ScaledStim, build_fiber  # noqa: E402

for diam in [3]:
    if diam == 3:
        nsim = 0
    else:
        nsim = 1
    for typ in ['twood', 'threed']:
        if typ == 'twood':
            continue  # temp
            potfile = rf'D:\threed_ascent\samples\{sample}2\models\0\sims\{sim}\n_sims\{nsim}\data\inputs\inner{innernum}_fiber{fibernum}.dat'
            pots = np.loadtxt(potfile, skiprows=1) * 1000  # mV
        else:
            potfile = rf'D:\threed_ascent\samples\{sample}3\models\0\sims\{sim}\n_sims\{nsim}\data\inputs\inner0_fiber{mfi}.dat'
            pots = np.loadtxt(potfile, skiprows=1)[::-1] * 1000  # mV

        n_sections = len(pots)

        model = FiberModel.MRG_INTERPOLATION

        # create fiber
        fiber = build_fiber(diameter=diam, fiber_model=model, temperature=37, n_sections=n_sections)
        fiber.potentials = pots

        # create biphasic square wave
        waveform = np.concatenate((np.ones(250), -np.ones(250), np.zeros(49600)))

        # parameters
        time_step = 0.001
        time_stop = 10

        # Create instance of ScaledStim class
        stimulation = ScaledStim(waveform=waveform, dt=time_step, tstop=time_stop)

        fiber.set_save_gating()
        fiber.set_save_vm()

        # run threshold search
        amp, ap = stimulation.find_threshold(fiber, stimamp_top=-5, stimamp_bottom=-0.1)

        plt.figure()
        plt.plot([apc.time for apc in fiber.apc])
        plt.title(f'Ap initation times for {typ} {diam} micron fiber')
        # %% make video of vm
        fps = 30
        skip = 5

        duration = len(fiber.vm[1]) / fps / skip
        ylim = (np.amin([v for v in fiber.vm[1:-1]]), np.amax([v for v in fiber.vm[1:-1]]))

        fig, ax = plt.subplots()

        def make_frame(i):
            ind = int(i * skip * fps)
            ax.clear()
            ax.set_ylim(ylim)
            ax.plot(fiber.coordinates[11:-11:11], [v[ind] for v in fiber.vm[1:-1]], lw=3)
            plt.title(f'Time: {i*skip / fps} ms')
            plt.ylim([-100, -40])
            return mplfig_to_npimage(fig)

        animation = VideoClip(make_frame, duration=duration)
        animation.write_videofile(f'{typ}_{diam}_sim{sim}.mp4', fps=fps)

        duration = len(fiber.vm[1]) / fps / skip
        ylim = (0, 1)
        #%%
        plt.figure()
        plt.plot(fiber.coordinates[::11][1:-1], -np.diff(fiber.potentials[::11], n=2))
        plt.gca().twinx().plot(fiber.coordinates[::11], [apc.time for apc in fiber.apc], 'r')
        test = [apc.time for apc in fiber.apc]
        test = np.array(test)
        test[test == 0] = float('Inf')
        print(np.argmin(test))
        plt.axvline(fiber.coordinates[::11][np.argmin(test)], color='k')

        apctest = np.loadtxt(
            r'D:\threed_ascent\samples\673\models\0\sims\3\n_sims\0\data\outputs\ap_loctime_inner0_fiber212_amp0.dat'
        )
        apctest = np.array(apctest)
        apctest[apctest == 0] = float('Inf')
        print(np.argmin(apctest))
        plt.xlim(15000, 35000)
        plt.axvline(fiber.coordinates[::11][np.argmin(apctest)], color='g')

        #%% plot ve over time, with each line colored according to time, downsample time
        plt.figure()
        vm = np.array([np.array(v) for v in fiber.vm[1:-1]]).T
        start = 1
        stop = 800
        step = 75
        times = np.arange(start, stop, step)
        for i, t in enumerate(times):
            plt.plot(fiber.coordinates[11:-11:11], vm[t], color=plt.cm.viridis(i / len(times)))
        plt.title(f'Vm over time for {typ} {diam} micron fiber')
        # make colorbar
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=start, vmax=stop))
        sm._A = []
        plt.colorbar(sm)
        plt.ylim(-100, -40)
        #%%
        ##% that didnt look good, let's try the same as above but make subplots, where each time is its own subplot
        fig, axs = plt.subplots(len(times), 1, sharex=True, sharey=True)
        for i, t in enumerate(times):
            axs[i].set_ylabel(f'{round(stimulation.time[t],2)} ms\nVm', color='r')
            # add second difference as black line on twinx behind the plot
            axs[i].twinx().plot(
                fiber.coordinates[11:-11:11],
                -np.diff(fiber.potentials[::11] * stimulation.waveform[t], n=2),
                color='k',
                linewidth=2,
            )
            axs[i].plot(fiber.coordinates[11:-11:11], vm[t], color='r', linewidth=2)
            axs[i].set_ylim(-90, -50)
            axs[i].twinx().set_ylabel('$d^2V_e$')

        # fig.suptitle(f'Vm over time for {typ} {diam} micron fiber')
        # make colorbar
        # sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=start, vmax=stop))
        # sm._A = []
        # plt.colorbar(sm)
        fig.set_size_inches(10, 20)
        plt.subplots_adjust(hspace=0)
