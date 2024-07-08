"""Defines main for managing NEURON simulations.

The copyrights of this software are owned by Duke University. Please
refer to the LICENSE and README.md files for licensing instructions. The
source code can be found on the following GitHub repository:
https://github.com/wmglab-duke/ascent
"""

import json
import sys
import time

from pyfibers import FiberModel, ScaledStim, _Fiber, build_fiber
from pyfibers.enums import BoundsSearchMode, TerminationMode, ThresholdCondition
from saving import Saving


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
    bounds_search_mode = getattr(BoundsSearchMode, bounds_search_configs['mode'])
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
    with open(waveform_path) as waveform_file:
        dt = float(waveform_file.readline().strip())  # time step
        tstop = int(waveform_file.readline().strip())  # stop time
        file_lines = waveform_file.read().splitlines()
        waveform = [float(i) for i in file_lines]

    # Read in extracellular potentials
    with open(potentials_path) as potentials_file:
        axontotal = int(potentials_file.readline())
        file_lines = potentials_file.read().splitlines()
        potentials = [float(i) * 1000 for i in file_lines]  # Need to convert to V -> mV

    # create fiber object
    model = getattr(FiberModel, sim_configs['fibers']['mode'])
    diameter = sim_configs['fibers']['z_parameters']['diameter']
    fiber = build_fiber(diameter=diameter, fiber_model=model, temperature=temperature, n_sections=axontotal)
    fiber.potentials = potentials

    # create stimulation object
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

    # create saving object
    saving_configs = sim_configs['saving']
    saving = handle_saving(fiber, fiber_ind, inner_ind, saving_configs, sim_path, stimulation.dt)

    # determine protocol
    protocol_configs = sim_configs['protocol']
    amps = protocol_configs['amplitudes'] if protocol_configs['mode'] == 'FINITE_AMPLITUDES' else False

    if not amps:  # enter binary search modes
        amp = threshold_protocol(fiber, protocol_configs, sim_configs, stimulation)
        saving.save_thresh(amp)  # Save threshold value to file
        saving.save_variables(fiber, stimulation)  # Save user-specified variables
        saving.save_runtime(time.time() - start_time)  # Save runtime of simulation

    else:  # finite amplitudes protocol
        time_total = 0
        for amp_ind, amp in enumerate(amps):
            print(f'Running amp {amp_ind} of {len(amps)}: {amp} mA')

            n_aps = stimulation.run_sim(stimamp=amp, fiber=fiber)
            time_individual = time.time() - start_time - time_total
            saving.save_variables(fiber, stimulation, amp_ind)  # Save user-specified variables
            saving.save_activation(n_aps, amp_ind)  # Save number of APs triggered
            saving.save_runtime(time_individual, amp_ind)  # Save runtime of inidividual run

            time_total += time_individual


def threshold_protocol(fiber, protocol_configs: dict, sim_configs: dict, stimulation: ScaledStim):
    """Prepare for bisection threshold search.

    :param fiber: fiber object
    :param protocol_configs: dictionary containing protocol configs from <sim_index>.json
    :param sim_configs: dictionary containing simulation configs from <sim_index>.json
    :param stimulation: instance of ScaledStim class
    :return: returns threshold amplitude (nA)
    """
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
    )
    print(f'Threshold found! {amp} mA for a fiber with diameter {sim_configs["fibers"]["z_parameters"]["diameter"]}')
    return amp


def handle_saving(
    fiber: _Fiber, fiber_ind: int, inner_ind: int, saving_configs: dict, sim_path: str, time_step: float
) -> Saving:
    """Create an instance of the Saving class.

    :param fiber: instance of Fiber class
    :param fiber_ind: index of fiber #
    :param inner_ind: index of inner #
    :param saving_configs: dictionary containing saving configs from <sim_index>.json
    :param sim_path: path to n_sim directory
    :param time_step: time step for simulation
    :return: returns instance of the Saving class
    """
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
        time_step,
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
    return saving


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
