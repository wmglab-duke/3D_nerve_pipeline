"""Defines main for managing NEURON simulations.

The copyrights of this software are owned by Duke University. Please
refer to the LICENSE and README.md files for licensing instructions. The
source code can be found on the following GitHub repository:
https://github.com/wmglab-duke/asc ent
"""

import json
import sys
import time

from saving import Saving
from wmglab_neuron import FiberModel, ScaledStim, build_fiber
from wmglab_neuron.enums import BoundsSearchMode, TerminationMode, ThresholdCondition


def handle_termination(protocol_configs: dict) -> (int, float):
    """Handle termination mode and tolerance.

    :param protocol_configs: dictionary containing protocol configs
    :return: returns termination mode and tolerance required for NEURON simulation
    """
    if protocol_configs['termination_criteria']['mode'] == 'PERCENT_DIFFERENCE':
        termination_mode = TerminationMode.PERCENT_DIFFERENCE
        termination_tolerance = protocol_configs['termination_criteria']['percent']
    elif protocol_configs['termination_criteria']['mode'] == 'ABSOLUTE_DIFFERENCE':
        termination_mode = TerminationMode.ABSOLUTE_DIFFERENCE
        termination_tolerance = protocol_configs['termination_criteria']['tolerance']
    else:
        return None
    return termination_mode, termination_tolerance


def handle_bounds_search(bounds_search_configs: dict) -> (str, float, float, float, int):
    """Handle bounds search configs for simulation.

    :param bounds_search_configs: dictionary containing information about bounds search for simulation
    :return: returns the values required for bounds searching during NEURON simulation
    """
    bounds_search_mode = bounds_search_configs['mode']
    bounds_search_step = bounds_search_configs['step']
    stimamp_top = bounds_search_configs['top']
    stimamp_bottom = bounds_search_configs['bottom']
    max_iterations = bounds_search_configs.get("max_steps", 100)
    return bounds_search_mode, bounds_search_step, stimamp_top, stimamp_bottom, max_iterations


def main(
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

    # Read in <sim_index>.json file as dictionary
    with open(f'{sim_path}/{n_sim}.json') as file:
        sim_configs = json.load(file)

    # Read in <model_index>.json file to get temperature
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

    # Determine fiber model
    fiber_model = sim_configs['fibers']['mode']
    model = getattr(FiberModel, fiber_model)

    # create fiber
    fiber = build_fiber(
        diameter=sim_configs['fibers']['z_parameters']['diameter'],
        fiber_model=model,
        temperature=float(temperature),
        n_sections=axontotal,
    )

    fiber.potentials = potentials

    # create stimulation
    stimulation = ScaledStim(waveform=waveform, dt=dt, tstop=tstop)

    # attach intracellular stimulation
    istim_configs = sim_configs['intracellular_stim']
    stimulation.set_intracellular_stim(
        delay=istim_configs['times']['IntraStim_PulseTrain_delay'],
        pw=istim_configs['times']['pw'],
        dur=istim_configs['times']['IntraStim_PulseTrain_dur'],
        freq=istim_configs['pulse_repetition_freq'],
        amp=istim_configs['amp'],
        ind=istim_configs['ind'],
    )

    # Determine saving protocols
    saving_configs = sim_configs['saving']

    # Determine optional saving configurations
    if saving_configs['space']['vm'] or saving_configs['time']['vm']:
        fiber.set_save_vm()
    if saving_configs['space']['gating'] or saving_configs['time']['gating']:
        fiber.set_save_gating()

    end_ap_times = 'end_ap_times' in saving_configs  # end_ap_times in <sim_index>.json
    loc_min = saving_configs['end_ap_times']['loc_min'] if end_ap_times else None
    loc_max = saving_configs['end_ap_times']['loc_max'] if end_ap_times else None
    ap_loctime = bool('aploctime' in saving_configs and saving_configs['aploctime'])  # ap_loc_time in <sim_ind>.json
    runtimes = bool('runtimes' in saving_configs and saving_configs['runtimes'])  # save runtimes

    # create Saving instance
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
        ap_loctime=ap_loctime,
        runtime=runtimes,
    )

    # determine protocol
    protocol_configs = sim_configs['protocol']
    amps = protocol_configs['amplitudes'] if protocol_configs['mode'] == 'FINITE_AMPLITUDES' else False

    if not amps:  # enter binary search modes
        ap_detect_location = protocol_configs['threshold']['ap_detect_location']
        istim_delay = istim_configs['times']['IntraStim_PulseTrain_delay']

        if protocol_configs['mode'] == 'ACTIVATION_THRESHOLD':
            condition = ThresholdCondition.ACTIVATION
        elif protocol_configs['mode'] == 'BLOCK_THRESHOLD':
            condition = ThresholdCondition.BLOCK

        # determine termination protocols for binary search
        termination_mode, termination_tolerance = handle_termination(protocol_configs)
        if 'bounds_search' not in protocol_configs:
            bounds_search_mode, bounds_search_step, stimamp_top, stimamp_bottom, max_iterations = (
                BoundsSearchMode.PERCENT_INCREMENT,
                10,
                -1,
                -0.01,
                100,
            )
        else:
            bounds_search_mode, bounds_search_step, stimamp_top, stimamp_bottom, max_iterations = handle_bounds_search(
                protocol_configs['bounds_search']
            )

        exit_t_shift = protocol_configs.get('exit_t_shift', 5)

        kwargs = {
            "ap_detect_location": ap_detect_location,
            "istim_delay": istim_delay,
        }

        # submit fiber for simulation
        amp, ap = stimulation.find_threshold(
            fiber,
            condition=condition,
            bounds_search_mode=bounds_search_mode,
            bounds_search_step=bounds_search_step,
            termination_mode=termination_mode,
            termination_tolerance=termination_tolerance,
            stimamp_top=stimamp_top,
            stimamp_bottom=stimamp_bottom,
            max_iterations=max_iterations,
            exit_t_shift=exit_t_shift,
            **kwargs,
        )
        print(f'Threshold found! {amp}nA for a fiber with diameter {sim_configs["fibers"]["z_parameters"]["diameter"]}')
        saving.save_thresh(amp)  # Save threshold value to file
        time_individual = time.time() - start_time
        saving.save_variables(fiber, stimulation)  # Save user-specified variables
        saving.save_runtime(time_individual)  # Save runtime of simulation

    else:  # finite amplitudes protocol
        time_total = 0
        for amp_ind, amp in enumerate(amps):
            print(f'Running amp {amp_ind} of {len(amps)}: {amp} mA')

            n_aps = stimulation.run_sim(amp, fiber)
            time_individual = time.time() - start_time - time_total
            saving.save_variables(fiber, stimulation, amp_ind)  # Save user-specified variables
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
