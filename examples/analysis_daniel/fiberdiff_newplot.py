import sys

import numpy as np
from shapely.geometry import Point

sys.path.append(r'C:\nrn\lib\python')  # noqa: E800
import os

os.chdir('../..')

import matplotlib.pyplot as plt
from moviepy.video.io.bindings import mplfig_to_npimage
from moviepy.video.VideoClip import VideoClip

sim = 3

diams = [3, 13]

n_sims = [0, 5]

sample = 57

innernum, fibernum, mfi = 11, 1, 26

samplename = '5R'

# %%
# sys.exit('temp')
# now simulate the 3 micron fiber 2D
from wmglab_neuron import FiberModel, ScaledStim, build_fiber  # noqa: E402

for diam in [3, 13]:
    if diam == 3:
        nsim = 0
    else:
        nsim = 5
    for typ in ['twood', 'threed']:
        if typ == 'twood':
            potfile = rf'D:\threed_ascent\samples\{sample}2\models\0\sims\{sim}\n_sims\{nsim}\data\inputs\inner{innernum}_fiber{fibernum}.dat'
        else:
            potfile = rf'D:\threed_ascent\samples\{sample}3\models\0\sims\{sim}\n_sims\{nsim}\data\inputs\inner0_fiber{mfi}.dat'

        pots = np.loadtxt(potfile, skiprows=1) * 1000  # mV

        n_sections = len(pots)

        model = FiberModel.MRG_INTERPOLATION

        # create fiber
        fiber = build_fiber(diameter=diam, fiber_model=model, temperature=37, n_sections=n_sections)
        fiber.potentials = pots

        # create biphasic square wave
        waveform = np.concatenate((np.ones(250), np.zeros(49600)))

        # parameters
        time_step = 0.001
        time_stop = 5

        # Create instance of ScaledStim class
        stimulation = ScaledStim(waveform=waveform, dt=time_step, tstop=time_stop)

        fiber.set_save_gating()
        fiber.set_save_vm()

        # run threshold search
        amp, ap = stimulation.find_threshold(fiber, stimamp_top=-5, stimamp_bottom=-0.1)

        # %% plot ve over time, with each line colored according to time, downsample time
        plt.figure()
        vm = np.array([np.array(v) for v in fiber.vm[1:-1]]).T
        start = 1
        stop = 400
        step = 25
        times = [5, 20, 50, 100, 300, 600]
        for i, t in enumerate(times):
            plt.plot(fiber.coordinates[11:-11:11], vm[t], color=plt.cm.viridis(i / len(times)))
        plt.title(f'Vm over time for 3 micron fiber')
        # make colorbar
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=start, vmax=stop))
        sm._A = []
        plt.colorbar(sm)
        plt.ylim(-100, -40)
        # %%
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
        plt.suptitle(f'{typ}, {diam}')
