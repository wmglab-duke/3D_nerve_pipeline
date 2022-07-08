import pickle
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '../')))
from src.utils import (Config, Configurable, DiamDistMode, Exceptionable, FiberGeometry,
                       FiberXYMode, FiberZMode, MyelinatedSamplingType, MyelinationMode, Saveable,
                       SetupMode, WriteMode, NeuronRunMode, TerminationCriteriaMode, SearchAmplitudeIncrementMode)
from stimulation import Stimulation
from saving import Saving

def main(fiber_path, inner_ind, fiber_ind, potentials_path, waveform_path, sim_path):
    # load in fiber object
    fiber = pickle.load(open(fiber_path, 'rb'))

    sim_config_path = os.path.join(sim_path, '0.json')

    # create stimulation object instance
    stim = Stimulation()
    stim \
        .add(SetupMode.NEW, Config.SIM, sim_config_path) \
        .load_potentials(potentials_path) \
        .load_waveform(waveform_path)

    # create necessary variables for simulation
    potentials = stim.potentials
    waveform = stim.waveform
    dt = stim.dt
    tstop = stim.tstop
    n_fiber_coords = len(stim.potentials)
    n_tsteps = len(waveform)

    # create NEURON fiber sections
    fiber.generate(n_fiber_coords)

    # attach intracellular stimulation
    if fiber.myelination:
        stim.apply_intracellular(fiber.node)
    else:
        stim.apply_intracellular(fiber.sec)

    # create saving object instance
    saving = Saving()
    saving \
        .add(SetupMode.NEW, Config.SIM, sim_config_path) \
        .inherit(sim_path, fiber.axonnodes, fiber.delta_z)

    # submit fiber for simulation
    fiber.submit(potentials, waveform, n_tsteps, dt, tstop, inner_ind, fiber_ind, saving)

    # save data
    saving.write2file(inner_ind, fiber_ind)

# load in arguments from command line
if __name__ == "__main__":  # Allows for the safe importing of the main module
    fiber_path = sys.argv[1]
    inner_ind = sys.argv[2]
    fiber_ind = sys.argv[3]
    potentials_path = sys.argv[4]
    waveform_path = sys.argv[5]
    sim_path = sys.argv[6]

    main(fiber_path, inner_ind, fiber_ind, potentials_path, waveform_path, sim_path)
    print('done with local_run_controls.py')



