"""Defines main for managing NEURON simulations.

The copyrights of this software are owned by Duke University. Please
refer to the LICENSE and README.md files for licensing instructions. The
source code can be found on the following GitHub repository:
https://github.com/wmglab-duke/ascent
"""

import pickle
import sys
import time
import json

from wmglab_neuron import Recording, Saving, Stimulation


def main(
        fiber_path: str,
        inner_ind: int,
        fiber_ind: int,
        potentials_path: str,
        waveform_path: str,
        sim_path: str,
        n_sim: int):
    """Control flow of a single n_sim NEURON simulation.

    :param fiber_path: path to fiber.obj
    :param inner_ind: inner index
    :param fiber_ind: fiber index
    :param potentials_path: path to potentials file
    :param waveform_path: path to waveform file
    :param sim_path: path to n_sim directory
    :param n_sim: n_sim number
    """
    start_time = time.time()  # Starting time of simulation

    # load in fiber object
    with open(fiber_path, 'rb') as fiber_file:
        fiber = pickle.load(fiber_file)
    fiber.inner_ind, fiber.fiber_ind = inner_ind, fiber_ind

    with open(f'{sim_path}/{n_sim}.json') as file:
        sim_configs = json.load(file)

    protocol_configs = sim_configs['protocol']
    saving_configs = sim_configs['saving']
    istim_configs = sim_configs['intracellular_stim']

    # Read in waveform array, time step, and stop time
    with open(waveform_path, 'r') as waveform_file:
        dt = float(waveform_file.readline().strip())  # time step
        tstop = int(waveform_file.readline().strip())  # stop time
        file_lines = waveform_file.read().splitlines()
        waveform = [float(i) for i in file_lines]

    # Read in extracellular potentials
    with open(potentials_path, 'r') as potentials_file:
        axontotal = int(potentials_file.readline())
        file_lines = potentials_file.read().splitlines()
        potentials = [float(i) * 1000 for i in file_lines]  # Need to convert to V -> mV

    stimulation = Stimulation(potentials_list=potentials, waveform_list=waveform, dt=dt, tstop=tstop)

    # create NEURON fiber sections
    fiber.generate(axontotal)

    # attach intracellular stimulation
    stimulation.apply_intracellular(fiber,
                                    delay=istim_configs['times']['IntraStim_PulseTrain_delay'],
                                    pw=istim_configs['times']['pw'],
                                    dur=istim_configs['times']['IntraStim_PulseTrain_dur'],
                                    freq=istim_configs['pulse_repetition_freq'],
                                    amp=istim_configs['amp'],
                                    ind=istim_configs['ind'],
                                    )

    # Determine optional saving configurations
    end_ap_times = True if 'end_ap_times' in saving_configs.keys() else False
    loc_min = saving_configs['end_ap_times']['loc_min'] if end_ap_times else None
    loc_max = saving_configs['end_ap_times']['loc_max'] if end_ap_times else None
    ap_end_thresh = saving_configs['end_ap_times']['threshold'] if end_ap_times else None
    ap_loctime = True if 'aploctime' in saving_configs.keys() else False
    runtime = True if 'runtime' in saving_configs.keys() else False

    # create saving object instance
    saving = Saving(stimulation.dt,
                    fiber,
                    space_vm=saving_configs['space']['vm'],
                    space_gating=saving_configs['space']['gating'],
                    space_times=saving_configs['space']['times'],
                    time_vm=saving_configs['time']['vm'],
                    time_gating=saving_configs['time']['gating'],
                    istim=saving_configs['time']['istim'],
                    locs=saving_configs['time']['locs'],
                    end_ap_times=end_ap_times,
                    loc_min=loc_min,
                    loc_max=loc_max,
                    ap_end_thresh=ap_end_thresh,
                    ap_loctime=ap_loctime,
                    runtime=runtime,
                    output_path=os.path.join(sim_path, 'data', 'outputs')
                    )

    # create recording object instance
    recording = Recording(fiber)

    # submit fiber for simulation
    amps = protocol_configs['amplitudes'] if protocol_configs['mode'] == 'FINITE_AMPLITUDES' else False
    fiber.submit(stimulation,
                 saving,
                 recording,
                 start_time,
                 protocol_mode=protocol_configs['mode'],
                 amps=amps,
                 t_init_ss=protocol_configs['initSS'],
                 dt_init_ss=protocol_configs['dt_initSS'],
                 )


# load in arguments from command line
if __name__ == "__main__":  # Allows for the safe importing of the main module
    fiber_path = sys.argv[1]
    inner_ind = sys.argv[2]
    fiber_ind = sys.argv[3]
    potentials_path = sys.argv[4]
    waveform_path = sys.argv[5]
    sim_path = sys.argv[6]
    n_sim = sys.argv[7]

    main(fiber_path, inner_ind, fiber_ind, potentials_path, waveform_path, sim_path, n_sim)
    print('done')
