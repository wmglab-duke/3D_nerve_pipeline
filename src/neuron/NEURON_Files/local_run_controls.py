import pickle
import sys
import time

from stimulation import Stimulation
from saving import Saving
from recording import Recording

def main(fiber_path, inner_ind, fiber_ind, potentials_path, waveform_path, sim_path):
    start_time = time.time()

    # load in fiber object
    fiber = pickle.load(open(fiber_path, 'rb'))
    fiber.inner_ind, fiber.fiber_ind = inner_ind, fiber_ind

    # create stimulation object instance
    stimulation = Stimulation()
    stimulation \
        .load_potentials(potentials_path) \
        .load_waveform(waveform_path)

    # create NEURON fiber sections
    n_fiber_coords = len(stimulation.potentials)
    fiber.generate(n_fiber_coords)

    # attach intracellular stimulation
    stimulation.apply_intracellular(fiber)

    # create saving object instance
    saving = Saving()
    saving.inherit(sim_path, stimulation.dt, fiber)

    # create recording object instance
    recording = Recording(fiber)

    # submit fiber for simulation
    fiber.submit(stimulation, saving, recording)

    # save data
    runtime = time.time()-start_time
    saving.write2file(recording, runtime, fiber, stimulation.dt)

# load in arguments from command line
if __name__ == "__main__":  # Allows for the safe importing of the main module
    fiber_path = sys.argv[1]
    inner_ind = sys.argv[2]
    fiber_ind = sys.argv[3]
    potentials_path = sys.argv[4]
    waveform_path = sys.argv[5]
    sim_path = sys.argv[6]

    main(fiber_path, inner_ind, fiber_ind, potentials_path, waveform_path, sim_path)



