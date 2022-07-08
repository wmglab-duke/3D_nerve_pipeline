import pickle
import os

from submit.python_stuff.stimulation import Stimulation

class RunControls():


fiber_path = os.path.join(os.getcwd(), 'n_sims', '0_0_0_0', 'fiber.obj')
fiber = pickle.load(open(fiber_path, 'rb'))

potentials_path = 'n_sims/0_0_0_0/data/inputs/inner0_fiber0.dat'
waveform_path = 'n_sims/0_0_0_0/data/inputs/waveform.dat'
stim = Stimulation()
stim.load_potentials(potentials_path)
stim.load_waveform(waveform_path)

potentials = stim.potentials
waveform = stim.waveform
dt = stim.dt
tstop = stim.tstop
n_fiber_coords = len(stim.potentials)
n_tsteps = len(waveform)

fiber.generate(n_fiber_coords)
fiber.submit(potentials, waveform, n_tsteps, dt, tstop)