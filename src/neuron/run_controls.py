"""Defines main for managing NEURON simulations.

The copyrights of this software are owned by Duke University. Please
refer to the LICENSE and README.md files for licensing instructions. The
source code can be found on the following GitHub repository:
https://github.com/wmglab-duke/ascent
"""

import json
import sys
import time

import numpy as np
from pyfibers import BoundsSearchMode, FiberModel, ScaledStim, TerminationMode, ThresholdCondition, build_fiber
from saving import initialize_saving, save_activation, save_runtime, save_thresh, save_variables


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


def handle_bounds_search(bounds_search_configs: dict) -> (str, float, int):
    """Handle bounds search configs for simulation.

    :param bounds_search_configs: dictionary containing information about bounds search for simulation
    :return: returns the values required for bounds searching during NEURON simulation
    """
    bounds_search_mode = getattr(BoundsSearchMode, bounds_search_configs['mode'])
    bounds_search_step = bounds_search_configs['step']
    max_iterations = bounds_search_configs.get("max_steps", 100)
    return bounds_search_mode, bounds_search_step, max_iterations


def main(
    inner_ind: int,
    fiber_ind: int,
    stimamp_top: float,
    stimamp_bottom: float,
    potentials_path: str,
    waveform_path: str,
    sim_path: str,
    n_sim: int,
):
    """Control flow of a single n_sim NEURON simulation.

    :param inner_ind: inner index
    :param fiber_ind: fiber index
    :param stimamp_top: top stimamp to start bounds search
    :param stimamp_bottom: bottom stimamp to start bounds search
    :param potentials_path: path to potentials file
    :param waveform_path: path to waveform file
    :param sim_path: path to n_sim directory
    :param n_sim: n_sim number
    :raises NotImplementedError: for deprecated features
    """
    start_time = time.time()  # Starting time of simulation
    stimamp_top, stimamp_bottom = float(stimamp_top), float(stimamp_bottom)

    # Read in <sim_index>.json file as dictionary
    with open(f'{sim_path}/{n_sim}.json') as file:
        sim_configs = json.load(file)

    # Read in <model_index>.json file to get temperature
    with open(f'{sim_path}/model.json') as file:
        model_configs = json.load(file)
        temperature = model_configs['temperature']

    # Read in waveform array, time step, and stop time
    with open(waveform_path) as waveform_file:
        dt = float(waveform_file.readline().strip())  # time step
        tstop = int(waveform_file.readline().strip())  # stop time
        file_lines = waveform_file.read().splitlines()
        waveform = np.array([float(i) for i in file_lines])

    # Read in extracellular potentials
    with open(potentials_path) as potentials_file:
        axontotal = int(potentials_file.readline())
        file_lines = potentials_file.read().splitlines()
        potentials = np.array([float(i) for i in file_lines]) * 1000  # Need to convert to V -> mV

    # create fiber object
    model = getattr(FiberModel, sim_configs['fibers']['mode'])
    diameter = sim_configs['fibers']['z_parameters']['diameter']
    fiber = build_fiber(diameter=diameter, fiber_model=model, temperature=temperature, n_sections=axontotal)
    fiber.potentials = potentials

    # attach intracellular stimulation
    # TODO change all references to istim in docs and code
    istim_configs = sim_configs.get('intrinsic_activity', {})
    if istim_configs:
        fiber.add_intrinsic_activity(**istim_configs)  # Add to docs to look at pyfiber docs for this

    if 'intracellular_stim' in istim_configs:
        raise NotImplementedError('Intracellular stimulation is deprecated, use intrinsic_activity instead')

    # initialize saving parameters
    saving_params = initialize_saving(sim_configs.get('saving', {}), fiber, fiber_ind, inner_ind, sim_path, dt)

    # determine protocol
    protocol_configs = sim_configs['protocol']
    t_init_ss, dt_init_ss = protocol_configs['initSS'], protocol_configs['dt_initSS']
    ap_detect_location = protocol_configs['threshold']['ap_detect_location']

    # create stimulation object
    stimulation = ScaledStim(waveform=waveform, dt=dt, tstop=tstop, t_init_ss=t_init_ss, dt_init_ss=dt_init_ss)

    if protocol_configs['mode'] != 'FINITE_AMPLITUDES':  # threshold search
        amp = threshold_protocol(
            fiber, protocol_configs, sim_configs, stimulation, stimamp_top, stimamp_bottom, ap_detect_location
        )
        save_thresh(saving_params, amp)  # Save threshold value to file
        save_variables(saving_params, fiber, stimulation)  # Save user-specified variables
        save_runtime(saving_params, time.time() - start_time)  # Save runtime of simulation

    else:  # finite amplitudes protocol
        time_total = 0
        amps = protocol_configs['amplitudes']
        for amp_ind, amp in enumerate(amps):
            print(f'Running amp {amp_ind} of {len(amps)}: {amp} mA')

            n_aps, _ = stimulation.run_sim(stimamp=amp, fiber=fiber, ap_detect_location=ap_detect_location)
            time_individual = time.time() - start_time - time_total
            save_variables(saving_params, fiber, stimulation, amp_ind)  # Save user-specified variables
            save_activation(saving_params, n_aps, amp_ind)  # Save number of APs triggered
            save_runtime(saving_params, time_individual, amp_ind)  # Save runtime of inidividual run

            time_total += time_individual


def threshold_protocol(
    fiber,
    protocol_configs: dict,
    sim_configs: dict,
    stimulation: ScaledStim,
    stimamp_top: float,
    stimamp_bottom: float,
    ap_detect_location: float,
) -> float:
    """Prepare for bisection threshold search.

    :param fiber: fiber object
    :param protocol_configs: dictionary containing protocol configs from <sim_index>.json
    :param sim_configs: dictionary containing simulation configs from <sim_index>.json
    :param stimulation: instance of ScaledStim class
    :param stimamp_top: top stimamp to start bounds search
    :param stimamp_bottom: bottom stimamp to start bounds search
    :param ap_detect_location: location to detect APs for threshold search
    :return: returns threshold amplitude (nA)
    """
    if protocol_configs['mode'] == 'ACTIVATION_THRESHOLD':
        condition = ThresholdCondition.ACTIVATION
    elif protocol_configs['mode'] == 'BLOCK_THRESHOLD':
        condition = ThresholdCondition.BLOCK
    # determine termination protocols for binary search
    termination_mode, termination_tolerance = handle_termination(protocol_configs)
    bounds_search_mode, bounds_search_step, max_iterations = handle_bounds_search(protocol_configs['bounds_search'])
    exit_t_shift = protocol_configs.get('exit_t_shift', 5)

    # submit fiber for simulation
    amp, _ = stimulation.find_threshold(
        fiber=fiber,
        condition=condition,
        bounds_search_mode=bounds_search_mode,
        bounds_search_step=bounds_search_step,
        termination_mode=termination_mode,
        termination_tolerance=termination_tolerance,
        stimamp_top=stimamp_top,
        stimamp_bottom=stimamp_bottom,
        max_iterations=max_iterations,
        exit_t_shift=exit_t_shift,
        ap_detect_location=ap_detect_location,
    )
    print(f'Threshold found! {amp} mA for a fiber with diameter {sim_configs["fibers"]["z_parameters"]["diameter"]}')
    return amp


# load in arguments from command line
if __name__ == "__main__":  # Allows for the safe importing of the main module
    main(*sys.argv[1:])
    print('done')
