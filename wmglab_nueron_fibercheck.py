"""Example use case of wmglab_neuron.

NOTE this is for development only
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage

sys.path.append(r'C:\nrn\lib\python')

from wmglab_neuron import FiberModel, ScaledStim, build_fiber  # noqa: E402

potfile = r'D:\threed_ascent\samples\250\models\0\sims\3\n_sims\5\data\inputs\inner0_fiber0.dat'

pots = np.loadtxt(potfile, skiprows=1) * 1000  # mV

n_sections = len(pots)

model = FiberModel.MRG_INTERPOLATION

# create fiber
fiber = build_fiber(diameter=13, fiber_model=model, temperature=37, n_sections=n_sections)
fiber.potentials = pots
# create curve of potentials
plt.plot(fiber.potentials)
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
amp, ap = stimulation.find_threshold(fiber, stimamp_top=-1, stimamp_bottom=-0.1)

plt.plot([apc.time for apc in fiber.apc])
#%% make video of vm
fps = 30
skip = 10

duration = len(fiber.vm[1]) / fps / skip
ylim = (np.amin([v for v in fiber.vm[1:-1]]), np.amax([v for v in fiber.vm[1:-1]]))

fig, ax = plt.subplots()


def make_frame(i):
    ind = int(i * skip * fps)
    ax.clear()
    ax.set_ylim(ylim)
    ax.plot(fiber.coordinates[11:-11:11], [v[ind] for v in fiber.vm[1:-1]], lw=3)
    plt.title(f'Time: {i/fps} ms')
    return mplfig_to_npimage(fig)


animation = VideoClip(make_frame, duration=duration)
animation.write_videofile('test13.mp4', fps=fps)

duration = len(fiber.vm[1]) / fps / skip
ylim = (0, 1)

#%%gating video
fig, ax = plt.subplots()


def make_frame(i):
    ind = int(i * skip * fps)
    ax.clear()
    ax.set_ylim(ylim)
    for var, data in fiber.gating.items():
        label = var if i == 0 else '_'
        ax.plot(fiber.coordinates[11:-11:11], [v[ind] for v in data[1:-1]], lw=2, label=var)
    plt.legend()
    plt.title(f'Time: {i/fps} ms')
    return mplfig_to_npimage(fig)


animation = VideoClip(make_frame, duration=duration)
animation.write_videofile('test13-gating.mp4', fps=fps)
#%% make video of vm at 110% threshold
stimulation.run_sim(amp * 1.1, fiber)

fps = 30
skip = 10

duration = len(fiber.vm[1]) / fps / skip
ylim = (np.amin([v for v in fiber.vm[1:-1]]), np.amax([v for v in fiber.vm[1:-1]]))

fig, ax = plt.subplots()


def make_frame(i):
    ind = int(i * skip * fps)
    ax.clear()
    ax.set_ylim(ylim)
    ax.plot(fiber.coordinates[11:-11:11], [v[ind] for v in fiber.vm[1:-1]], lw=3)
    plt.title(f'Time: {i/fps} ms')
    return mplfig_to_npimage(fig)


animation = VideoClip(make_frame, duration=duration)
animation.write_videofile('test13_1.1.mp4', fps=fps)
#%%now 3 micron

potfile = r'D:\threed_ascent\samples\250\models\0\sims\3\n_sims\0\data\inputs\inner0_fiber0.dat'

pots = np.loadtxt(potfile, skiprows=1) * 1000  # mV

n_sections = len(pots)

model = FiberModel.MRG_INTERPOLATION

# create fiber
fiber = build_fiber(diameter=3, fiber_model=model, temperature=37, n_sections=n_sections)
fiber.potentials = pots
# create curve of potentials
plt.plot(fiber.potentials)
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
amp, ap = stimulation.find_threshold(fiber, stimamp_top=-1, stimamp_bottom=-0.1)

plt.figure()
plt.plot([apc.time for apc in fiber.apc])

#%% make video of vm
fps = 30
skip = 10

duration = len(fiber.vm[1]) / fps / skip
ylim = (np.amin([v for v in fiber.vm[1:-1]]), np.amax([v for v in fiber.vm[1:-1]]))

fig, ax = plt.subplots()


def make_frame(i):
    ind = int(i * skip * fps)
    ax.clear()
    ax.set_ylim(ylim)
    ax.plot(fiber.coordinates[11:-11:11], [v[ind] for v in fiber.vm[1:-1]], lw=3)
    plt.title(f'Time: {i/fps} ms')
    return mplfig_to_npimage(fig)


animation = VideoClip(make_frame, duration=duration)
animation.write_videofile('test3.mp4', fps=fps)
#%% make video of vm at 110% threshold
stimulation.run_sim(amp * 0.9, fiber)

fps = 30
skip = 10

duration = len(fiber.vm[1]) / fps / skip
ylim = (np.amin([v for v in fiber.vm[1:-1]]), np.amax([v for v in fiber.vm[1:-1]]))

fig, ax = plt.subplots()


def make_frame(i):
    ind = int(i * skip * fps)
    ax.clear()
    ax.set_ylim(ylim)
    ax.plot(fiber.coordinates[11:-11:11], [v[ind] for v in fiber.vm[1:-1]], lw=3)
    plt.title(f'Time: {i/fps} ms')

    return mplfig_to_npimage(fig)


animation = VideoClip(make_frame, duration=duration)
animation.write_videofile('test3_0.9.mp4', fps=fps)

#%%gating video
fig, ax = plt.subplots()
ylim = (0, 1)


def make_frame(i):
    ind = int(i * skip * fps)
    ax.clear()
    ax.set_ylim(ylim)
    for var, data in fiber.gating.items():
        label = var if i == 0 else '_'
        ax.plot(fiber.coordinates[11:-11:11], [v[ind] for v in data[1:-1]], lw=2, label=var)
    plt.title(f'Time: {i/fps} ms')
    plt.legend()
    return mplfig_to_npimage(fig)


animation = VideoClip(make_frame, duration=duration)
animation.write_videofile('test3-gating.mp4', fps=fps)
