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
# %% gating plot
fig, axs = plt.subplots(2, 1, sharex=True)
apc = [apc.time for apc in fiber.apc]
apc[0] = float('Inf')
apc[-1] = float('Inf')
actloc = np.argmin(apc)
# plot gating vars at activation location
for var, data in fiber.gating.items():
    plt.plot(stimulation.time, data[actloc], label=var)
plt.gca().twinx().plot(list(stimulation.time)[1:], -stimulation.waveform, 'k--')
plt.legend()
plt.xlim(0, 2)
plt.xlabel('ms')
axs[0].plot(stimulation.time, fiber.vm[actloc])
fig.set_size_inches([6, 6])
axs[0].set_title('Activation location 13 um')

fig, axs = plt.subplots(2, 1, sharex=True)
actloc = np.abs(len(fiber.nodes) - actloc)
# plot gating vars at activation location
for var, data in fiber.gating.items():
    plt.plot(stimulation.time, data[actloc], label=var)
plt.gca().twinx().plot(list(stimulation.time)[1:], stimulation.waveform, 'k--')
plt.legend()
plt.xlim(0, 2)
plt.xlabel('ms')
axs[0].plot(stimulation.time, fiber.vm[actloc])
fig.set_size_inches([6, 6])

# %% plot ve over time, with each line colored according to time, downsample time
plt.figure()
vm = np.array([np.array(v) for v in fiber.vm[1:-1]]).T
start = 1
stop = 800
step = 75
times = np.arange(start, stop, step)
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

# %%now 3 micron

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

# %% gating plot
fig, axs = plt.subplots(2, 1, sharex=True)
apc = [apc.time for apc in fiber.apc]
apc[0] = float('Inf')
apc[-1] = float('Inf')
actloc = np.argmin(apc)
# plot gating vars at activation location
for var, data in fiber.gating.items():
    plt.plot(stimulation.time, data[actloc], label=var)
plt.gca().twinx().plot(list(stimulation.time)[1:], -stimulation.waveform, 'k--')
plt.legend()
plt.xlim(0, 2)
plt.xlabel('ms')
axs[0].plot(stimulation.time, fiber.vm[actloc])
fig.set_size_inches([6, 6])
axs[0].set_title('Activation location 3 um')

fig, axs = plt.subplots(2, 1, sharex=True)
actloc = np.abs(len(fiber.nodes) - actloc)
# plot gating vars at activation location
for var, data in fiber.gating.items():
    plt.plot(stimulation.time, data[actloc], label=var)
plt.gca().twinx().plot(list(stimulation.time)[1:], stimulation.waveform, 'k--')
plt.legend()
plt.xlim(0, 2)
plt.xlabel('ms')
axs[0].plot(stimulation.time, fiber.vm[actloc])
fig.set_size_inches([6, 6])

# %% plot ve over time, with each line colored according to time, downsample time
plt.figure()
vm = np.array([np.array(v) for v in fiber.vm[1:-1]]).T
start = 1
stop = 800
step = 75
times = np.arange(start, stop, step)
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
