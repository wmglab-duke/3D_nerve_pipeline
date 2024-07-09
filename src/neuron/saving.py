"""Defines functinos for data saving after NEURON simulations.

The copyrights of this software are owned by Duke University. Please
refer to the LICENSE and README.md files for licensing instructions. The
source code can be found on the following GitHub repository:
https://github.com/wmglab-duke/ascent
"""

import os

import pandas as pd
from pyfibers import Fiber, ScaledStim


def initialize_saving(
    saving_configs: dict, fiber: Fiber, fiber_ind: int, inner_ind: int, sim_path: str, dt: float
) -> dict:
    """Initialize saving parameters.

    :param saving_configs: dictionary containing saving configurations (from sim.json)
    :param fiber: instance of Fiber class
    :param fiber_ind: index of fiber in list of fibers
    :param inner_ind: index of inner in list of inners
    :param sim_path: path to simulation directory
    :param dt: time step of stimulation (ms)
    :return: dictionary containing saving parameters
    """
    # Determine optional saving configurations and set fiber saving properties
    # TODO make saving and all sub arguments optional
    if saving_configs['space']['vm'] or saving_configs['time']['vm']:
        fiber.set_save_vm()
    if saving_configs['space']['gating'] or saving_configs['time']['gating']:
        fiber.set_save_gating()

    space_times = saving_configs['space']['times']
    locs = saving_configs['time']['locs']
    nodecount = len(fiber.nodes) if fiber else 0

    time_inds = [int(t / dt) for t in (space_times or [])]
    time_inds.sort()
    node_inds = [int((nodecount - 1) * loc) for loc in (locs or [])] if locs != 'all' else list(range(0, nodecount))
    node_inds.sort()
    output_path = os.path.join(sim_path, 'data', 'outputs')
    os.makedirs(output_path, exist_ok=True)

    return {
        'inner_ind': inner_ind,
        'fiber_ind': fiber_ind,
        'space_vm': saving_configs['space']['vm'],
        'space_gating': saving_configs['space']['gating'],
        'time_inds': time_inds,
        'time_vm': saving_configs['time']['vm'],
        'time_gating': saving_configs['time']['gating'],
        'istim': saving_configs['time']['istim'],
        'node_inds': node_inds,
        'ap_loctime': saving_configs.get('aploctime', False),
        'runtime': saving_configs.get('runtimes', False),
        'output_path': output_path,
    }


def save_thresh(params: dict, thresh: float):
    """Save threshold from NEURON simulation to file.

    :param params: dictionary containing saving parameters
    :param thresh: activation threshold from NEURON simulation
    """
    thresh_path = os.path.join(
        params['output_path'], f'thresh_inner{params["inner_ind"]}_fiber{params["fiber_ind"]}.dat'
    )
    with open(thresh_path, 'w') as thresh_file:
        thresh_file.write(f"{thresh:.6f}")


def save_runtime(params: dict, runtime: float, amp_ind: int = 0):
    """Save NEURON simulation runtime to file.

    :param params: dictionary containing saving parameters
    :param runtime: runtime of NEURON simulation
    :param amp_ind: index of amplitude in list of amplitudes for finite_amplitude protocol
    """
    if params['runtime']:
        runtimes_path = os.path.join(
            params['output_path'], f'runtime_inner{params["inner_ind"]}_fiber{params["fiber_ind"]}_amp{amp_ind}.dat'
        )
        with open(runtimes_path, 'w') as runtime_file:
            runtime_file.write(f'{runtime:.3f}')


def save_activation(params: dict, n_aps: int, amp_ind: int):
    """Save the number of action potentials that occurred at the location specified in sim config to file.

    :param params: dictionary containing saving parameters
    :param n_aps: number of action potentials that occurred
    :param amp_ind: index of amplitude
    """
    output_file_path = os.path.join(
        params['output_path'], f'activation_inner{params["inner_ind"]}_fiber{params["fiber_ind"]}_amp{amp_ind}.dat'
    )
    with open(output_file_path, 'w') as activation_file:
        activation_file.write(f'{n_aps}')


def save_variables(params: dict, fiber: Fiber, stimulation: ScaledStim, amp_ind: int = 0):
    """Write user-specified variables to file.

    :param params: dictionary containing saving parameters
    :param fiber: instance of Fiber class
    :param stimulation: instance of ScaledStim class
    :param amp_ind: index of amplitude if protocol is FINITE_AMPLITUDES
    """
    params['time_inds'] = [t for t in params['time_inds'] if t < len(stimulation.time)]
    if fiber.vm is not None:
        target_len = max(len(vm) for vm in fiber.vm if vm is not None)
        vm_data = pd.DataFrame(list(vm) if vm is not None else [None] * target_len for vm in fiber.vm)
        vm_data = vm_data.fillna(0)
        save_data(params, amp_ind, fiber, vm_data, stimulation.dt, 'vm', 'space', stimulation, 'mV')
        save_data(params, amp_ind, fiber, vm_data, stimulation.dt, 'vm', 'time', stimulation, 'mV')
    if fiber.gating is not None:
        all_gating_data = {}
        for gating_parameter, gating_data in fiber.gating.items():
            target_len = max(len(g) for g in gating_data if g is not None)
            all_gating_data[gating_parameter] = pd.DataFrame(
                [list(g) if g is not None else [None] * target_len for g in gating_data]
            )
            all_gating_data[gating_parameter] = all_gating_data[gating_parameter].fillna(0)
        save_gating_data(params, all_gating_data, amp_ind, fiber, stimulation.dt, 'space', stimulation)
        save_gating_data(params, all_gating_data, amp_ind, fiber, stimulation.dt, 'time', stimulation)
    if params['istim']:
        istim_data = pd.DataFrame(list(stimulation.istim_record))
        save_data(params, amp_ind, fiber, istim_data, stimulation.dt, 'istim', 'time', stimulation, 'nA')
    if params['ap_loctime']:
        save_aploctime(params, amp_ind, fiber)


def save_data(
    params: dict,
    amp_ind: int,
    fiber: Fiber,
    data: pd.DataFrame,
    dt: float,
    var_type: str,
    save_type: str,
    stimulation: ScaledStim,
    units: str = None,
):
    """General method to save data to file.

    :param params: dictionary containing saving parameters
    :param amp_ind: index of the amplitude in the finite amplitude protocol
    :param fiber: instance of Fiber class
    :param data: DataFrame containing the data to be saved
    :param dt: time step of stimulation (ms)
    :param var_type: type of variable to be saved (vm, istim, gating parameters)
    :param save_type: function of variable to be saved (space or time)
    :param stimulation: instance of ScaledStim class
    :param units: units of variable (mV, nA)
    """
    if (
        (save_type == 'space' and params['space_vm'])
        or (save_type == 'time' and params['time_vm'])
        or (var_type == 'istim' and params['istim'])
    ):
        file_path = os.path.join(
            params['output_path'],
            f'{var_type}_{save_type}_inner{params["inner_ind"]}_fiber{params["fiber_ind"]}_amp{amp_ind}.dat',
        )
        if save_type == 'space':
            data = data[params['time_inds']]  # save data only at user-specified times
            data.insert(0, 'Node#', [*range(0, len(fiber.nodes))])
        elif save_type == 'time':
            data = data.T[params['node_inds']]  # save data only at user-specified locations
            data.insert(0, 'Time', stimulation.time)
        header = handle_header(params, save_type=save_type, var_type=var_type, dt=dt, units=units)
        data.to_csv(file_path, header=header, sep=' ', float_format='%.6f', index=False)


def save_gating_data(
    params: dict, all_gating_data: dict, amp_ind: int, fiber: Fiber, dt: float, save_type: str, stimulation: ScaledStim
):
    """General method to save gating data to file.

    :param params: dictionary containing saving parameters
    :param all_gating_data: dictionary containing gating data DataFrames
    :param amp_ind: index of amplitude in finite amplitudes protocol
    :param fiber: instance of Fiber class
    :param dt: time step of simulation (ms)
    :param save_type: function of variable to be saved (space or time)
    :param stimulation: instance of ScaledStim class
    """
    if (save_type == 'space' and params['space_gating']) or (save_type == 'time' and params['time_gating']):
        for gating_param, gating_data in all_gating_data.items():
            save_data(params, amp_ind, fiber, gating_data, dt, "gating_" + gating_param, save_type, stimulation, '')


def handle_header(params: dict, save_type: str, var_type: str, dt: float = None, units: str = None):
    """Create a header for text file.

    :param params: dictionary containing saving parameters
    :param save_type: function of variable to be saved, can be function of time ('time') or space ('space')
    :param var_type: type of variable to be saved (Vm, h, mp, m, s)
    :param units: units of variable (mV, nA)
    :param dt: float value for stimulation time step (ms)
    :return: list of column headers
    """
    if save_type == 'space':  # F(x) - function of space
        return space_header(params, var_type=var_type, dt=dt, units=units)
    elif save_type == 'time':  # F(t) - function of time
        return time_header(params, var_type=var_type, units=units)
    return []


def space_header(params: dict, var_type: str, dt: float, units: str = None) -> list[str]:
    """Create header for F(x) files.

    :param params: dictionary containing saving parameters
    :param var_type: type of variable to be saved (Vm, h, mp, m, s)
    :param dt: float value of stimulation time step (ms)
    :param units: units of variable (mV, nA)
    :return: returns list of strings to be used in file header
    """
    header = ['Node#']
    suffix = f'({units})' if units else ''
    for time in params['time_inds']:
        header.append(f'{var_type}_time{(float(time) * float(dt))}ms{suffix}')
    return header


def time_header(params: dict, var_type: str, units: str = None) -> list[str]:
    """Create header for F(t) files.

    :param params: dictionary containing saving parameters
    :param var_type: type of variable to be saved (Vm, h, mp, m, s)
    :param units: units of variable (mV, nA)
    :return: returns list of strings to be used as file header
    """
    header = ['Time(ms)']
    suffix = f'({units})' if units else ''
    if var_type != 'istim':
        for node in params['node_inds']:
            header.append(f'{var_type}_node{node}{suffix}')
    elif var_type == 'istim':
        header.append(f'Istim({units})')
    return header


def save_aploctime(params: dict, amp_ind: int, fiber: Fiber):
    """Save time when last NoR AP was detected for each node to file.

    :param params: dictionary containing saving parameters
    :param amp_ind: index of amplitude in finite amplitudes protocol
    :param fiber: instance of Fiber class
    """
    aploctime_path = os.path.join(
        params['output_path'], f'ap_loctime_inner{params["inner_ind"]}_fiber{params["fiber_ind"]}_amp{amp_ind}.dat'
    )
    with open(aploctime_path, 'w') as file:
        for node_ind, _ in enumerate(fiber.nodes):
            file.write(f"{fiber.apc[node_ind].time:.6f}\n")
