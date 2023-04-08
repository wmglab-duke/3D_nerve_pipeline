"""Defines main for managing NEURON simulations.

The copyrights of this software are owned by Duke University. Please
refer to the LICENSE and README.md files for licensing instructions. The
source code can be found on the following GitHub repository:
https://github.com/wmglab-duke/ascent
"""

import json
import sys
import time

from saving import Saving

from wmglab_neuron import FiberBuilder, FiberModel, Stimulation
from wmglab_neuron.enums import TerminationMode, ThresholdCondition


def main(  # noqa: C901
    inner_ind: int,
    fiber_ind: int,
    potentials_path: str,
    waveform_path: str,
    sim_path: str,
    n_sim: int,
):
    """Control flow of a single n_sim NEURON simulation.

    :param inner_ind: inner index
    :param fiber_ind: fiber index
    :param potentials_path: path to potentials file
    :param waveform_path: path to waveform file
    :param sim_path: path to n_sim directory
    :param n_sim: n_sim number
    """
    start_time = time.time()  # Starting time of simulation

    with open(f'{sim_path}/{n_sim}.json') as file:
        sim_configs = json.load(file)

    with open(f'{sim_path}/model.json') as file:
        model_configs = json.load(file)

    temperature = model_configs['temperature']

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

    fiber_model = sim_configs['fibers']['mode']
    if fiber_model == 'MRG_DISCRETE':
        model = FiberModel.MRG_DISCRETE
    elif fiber_model == 'MRG_INTERPOLATION':
        model = FiberModel.MRG_INTERPOLATION
    elif fiber_model == 'RATTAY':
        model = FiberModel.RATTAY
    elif fiber_model == 'TIGERHOLM':
        model = FiberModel.TIGERHOLM
    elif fiber_model == 'SUNDT':
        model = FiberModel.SUNDT
    elif fiber_model == 'SCHILD94':
        model = FiberModel.SCHILD94
    elif fiber_model == 'SCHILD97':
        model = FiberModel.SCHILD97
    else:
        exit()

    # create fiber
    fiber = FiberBuilder.generate(
        diameter=sim_configs['fibers']['z_parameters']['diameter'],
        fiber_model=model,
        temperature=float(temperature),
        n_fiber_coords=axontotal,
    )

    protocol_configs = sim_configs['protocol']
    saving_configs = sim_configs['saving']
    istim_configs = sim_configs['intracellular_stim']

    stimulation = Stimulation(fiber, waveform=waveform, potentials=potentials, dt=dt, tstop=tstop)

    if saving_configs['space']['vm'] or saving_configs['time']['vm']:
        fiber.set_save_vm()
    if saving_configs['space']['gating'] or saving_configs['time']['gating']:
        fiber.set_save_gating()

    # attach intracellular stimulation
    stimulation.add_intracellular_stim(
        delay=istim_configs['times']['IntraStim_PulseTrain_delay'],
        pw=istim_configs['times']['pw'],
        dur=istim_configs['times']['IntraStim_PulseTrain_dur'],
        freq=istim_configs['pulse_repetition_freq'],
        amp=istim_configs['amp'],
        ind=istim_configs['ind'],
    )

    # Determine optional saving configurations
    end_ap_times = 'end_ap_times' in saving_configs
    loc_min = saving_configs['end_ap_times']['loc_min'] if end_ap_times else None
    loc_max = saving_configs['end_ap_times']['loc_max'] if end_ap_times else None
    ap_end_thresh = saving_configs['end_ap_times']['threshold'] if end_ap_times else None
    ap_loctime = bool('aploctime' in saving_configs and saving_configs['aploctime'])
    runtimes = bool('runtimes' in saving_configs and saving_configs['runtimes'])

    # create saving object instance
    if saving_configs['space']['vm'] is not None or saving_configs['time']['vm'] is not None:
        fiber.set_save_vm()
    if saving_configs['space']['gating'] is not None or saving_configs['time']['gating'] is not None:
        fiber.set_save_gating()

    saving = Saving(
        inner_ind,
        fiber_ind,
        sim_path,
        stimulation.dt,
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
        runtime=runtimes,
    )

    # submit fiber for simulation
    amps = protocol_configs['amplitudes'] if protocol_configs['mode'] == 'FINITE_AMPLITUDES' else False

    ap_detect_location = protocol_configs['threshold']['ap_detect_location']
    istim_delay = istim_configs['times']['IntraStim_PulseTrain_delay']
    run_kwargs = {
        'ap_detect_location': ap_detect_location,
        'istim_delay': istim_delay,
    }
    kwargs = {k: v for k, v in run_kwargs.items() if v is not None}

    if not amps:
        if protocol_configs['mode'] == 'ACTIVATION_THRESHOLD':
            condition = ThresholdCondition.ACTIVATION
        elif protocol_configs['mode'] == 'BLOCK_THRESHOLD':
            condition = ThresholdCondition.BLOCK
        if protocol_configs['termination_criteria']['mode'] == 'PERCENT_DIFFERENCE':
            termination_mode = TerminationMode.PERCENT_DIFFERENCE
            termination_tolerance = protocol_configs['termination_criteria']['percent']
        elif protocol_configs['termination_criteria']['mode'] == 'ABSOLUTE_DIFFERENCE':
            termination_mode = TerminationMode.ABSOLUTE_DIFFERENCE
            termination_tolerance = protocol_configs['termination_criteria']['tolerance']

        if 'bounds_search' in protocol_configs:
            bounds_search_mode = protocol_configs['bounds_search']['mode']
            bounds_search_step = protocol_configs['bounds_search']['step']
            stimamp_top = protocol_configs['bounds_search']['top']
            stimamp_bottom = protocol_configs['bounds_search']['bottom']
            max_iterations = protocol_configs['bounds_search'].get("max_steps")
        else:
            bounds_search_mode = None
            bounds_search_step = None
            stimamp_top = None
            stimamp_bottom = None

        # todo: implement check_threshold_interval

        threshold_args = [
            condition,
            bounds_search_mode,
            bounds_search_step,
            termination_mode,
            termination_tolerance,
            stimamp_top,
            stimamp_bottom,
            max_iterations,
        ]

        args = [arg for arg in threshold_args if arg is not None]

        amp, ap = stimulation.find_threshold(
            *args,
            **kwargs,
        )
        print(f'Threshold found! {amp}nA for a fiber with diameter {sim_configs["fibers"]["z_parameters"]["diameter"]}')
        saving.save_thresh(amp)  # Save threshold value to file
        time_individual = time.time() - start_time
        saving.save_variables(fiber, stimulation)  # Save user-specified variables
        saving.save_runtime(time_individual)  # Save runtime of simulation

    else:
        time_total = 0
        for amp_ind, amp in enumerate(amps):
            print(f'Running amp {amp_ind} of {len(amps)}: {amp} mA')

            n_aps = stimulation.run_sim(
                stimamp=amp,
                **kwargs,
            )
            time_individual = time.time() - start_time - time_total
            saving.save_variables(fiber, stimulation)  # Save user-specified variables
            saving.save_activation(n_aps, amp_ind)  # Save number of APs triggered
            saving.save_runtime(time_individual, amp_ind)  # Save runtime of inidividual run

            time_total += time_individual


# load in arguments from command line
if __name__ == "__main__":  # Allows for the safe importing of the main module
    inner_ind = sys.argv[1]
    fiber_ind = sys.argv[2]
    potentials_path = sys.argv[3]
    waveform_path = sys.argv[4]
    sim_path = sys.argv[5]
    n_sim = sys.argv[6]

    main(
        inner_ind,
        fiber_ind,
        potentials_path,
        waveform_path,
        sim_path,
        n_sim,
    )
    print('done')
