"""The copyrights of this software are owned by Duke University.

Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

import pickle
import sys
import time

from NEURON_Files import Recording, Saving, Stimulation


def main(fiber_path: str, inner_ind: int, fiber_ind: int, potentials_path: str, waveform_path: str, sim_path: str):
    """Control flow of a single n_sim NEURON simulation.

    :param fiber_path: path to fiber.obj
    :param inner_ind: inner index
    :param fiber_ind: fiber index
    :param potentials_path: path to potentials file
    :param waveform_path: path to waveform file
    :param sim_path: path to n_sim directory
    """
    start_time = time.time()  # Starting time of simulation

    # load in fiber object
    with open(fiber_path, 'rb') as fiber_file:
        fiber = pickle.load(fiber_file)
    fiber.inner_ind, fiber.fiber_ind = inner_ind, fiber_ind

    # create stimulation object instance
    stimulation = Stimulation()
    stimulation.load_potentials(potentials_path).load_waveform(waveform_path)

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
    fiber.submit(stimulation, saving, recording, start_time)


# load in arguments from command line
if __name__ == "__main__":  # Allows for the safe importing of the main module
    fiber_path = sys.argv[1]
    inner_ind = sys.argv[2]
    fiber_ind = sys.argv[3]
    potentials_path = sys.argv[4]
    waveform_path = sys.argv[5]
    sim_path = sys.argv[6]

    main(fiber_path, inner_ind, fiber_ind, potentials_path, waveform_path, sim_path)
    print('done')
