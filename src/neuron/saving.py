"""Defines Saving class.

The copyrights of this software are owned by Duke University. Please
refer to the LICENSE and README.md files for licensing instructions. The
source code can be found on the following GitHub repository:
https://github.com/wmglab-duke/ascent
"""

import os

import pandas as pd
from wmglab_neuron import ScaledStim, _Fiber


class Saving:
    """Manage saving parameters to file for NEURON simulations."""

    def __init__(
        self,
        inner_ind: int,
        fiber_ind: int,
        sim_path: str,
        dt: float,
        fiber: object,
        space_vm: bool = False,
        space_gating: bool = False,
        space_times: list = None,
        time_vm: bool = False,
        time_gating: bool = False,
        istim: bool = False,
        locs: list = None,
        end_ap_times: bool = False,
        loc_min: float = 0.1,
        loc_max: float = 0.9,
        ap_loctime: bool = False,
        runtime: bool = False,
    ):
        """Initialize instance of Saving class.

        :param inner_ind: inner #
        :param fiber_ind: fiber #
        :param sim_path: path to n_sim directory
        :param dt: user-specified time step for simulation (ms)
        :param fiber: instance of Fiber class
        :param space_vm: save transmembrane at all sections at the time stamps defined in 'space_times'
        :param space_gating: save channel gating parameters at all sections at the time stamps defined in 'space_times'
        :param space_times: times in the simulation at which to save the values of state variables [ms]
        :param time_vm: save the transmembrane potential at each time step at the locations defined in 'locs'
        :param time_gating: save the channel gating parameters at each time step at the locations defined in 'locs'
        :param istim: save the applied intracellular stimulation at each time step
        :param locs: locations (decimal percentages of fiber length) at which to save state variables at all time steps
        :param end_ap_times: record when action potential occurs at specified indices
        :param loc_min: if end_ap_times, decimal % of fiber length at which to save times when AP is triggered
        :param loc_max: if end_ap_times, decimal % of fiber length at which to save times when Vm returns to threshold
        :param ap_loctime: save, for each fiber node, the last time an AP passed over that node
        :param runtime: save the simulation runtime
        :return: Saving object
        """
        nodecount = len(fiber.nodes)
        self.inner_ind = inner_ind
        self.fiber_ind = fiber_ind
        self.space_vm = space_vm
        self.space_gating = space_gating
        sim_times = space_times
        self.time_inds = [int(t / dt) for t in sim_times]  # divide by dt to avoid indexing error
        self.time_inds.sort()
        self.time_vm = time_vm
        self.time_gating = time_gating
        self.istim = istim
        self.locs = locs
        if self.locs != 'all':
            self.node_inds = [int((nodecount - 1) * loc) for loc in self.locs]
        elif self.locs == 'all':
            self.node_inds = list(range(0, nodecount))
        self.node_inds.sort()
        if end_ap_times:
            self.ap_end_times = True
            node_ind_min = int((nodecount - 1) * loc_min)
            node_ind_max = int((nodecount - 1) * loc_max)
            self.ap_end_inds = [node_ind_min, node_ind_max]
        else:
            self.ap_end_times = None
        self.ap_loctime = ap_loctime
        self.runtime = runtime
        self.output_path = os.path.join(sim_path, 'data', 'outputs')
        os.makedirs(self.output_path, exist_ok=True)
        return

    def save_thresh(self, thresh: float):
        """Save threshold from NEURON simulation to file.

        :param thresh: activation threshold from NEURON simulation
        """
        # save threshold to submit/n_sims/#/data/outputs/thresh_inner#_fiber#.dat
        thresh_path = os.path.join(self.output_path, f'thresh_inner{self.inner_ind}_fiber{self.fiber_ind}.dat')
        with open(thresh_path, 'w') as thresh_file:
            thresh_file.write(f"{thresh:.6f} mA")

    def save_runtime(self, runtime: float, amp_ind: int = 0):
        """Save NEURON simulation runtime to file.

        :param runtime: runtime of NEURON simulation
        :param amp_ind: index of amplitude in list of amplitudes for finite_amplitude protocol
        """
        if self.runtime:
            # save runtime to submit/n_sims/#/data/outputs/runtime_inner#_fiber#_amp#.dat
            runtimes_path = os.path.join(
                self.output_path, f'runtime_inner{self.inner_ind}_fiber{self.fiber_ind}_amp{amp_ind}.dat'
            )
            with open(runtimes_path, 'w') as runtime_file:
                runtime_file.write(f'{runtime:.3f}s')

    def save_activation(self, n_aps: int, amp_ind: int):
        """Save the number of action potentials that occurred at the location specified in sim config to file.

        :param n_aps: number of action potentials that occurred
        :param amp_ind: index of amplitude
        """
        # save number of action potentials to submit/n_sims/#/data/outputs/activation_inner#_fiber#_amp#.dat
        output_file_path = os.path.join(
            self.output_path, f'activation_inner{self.inner_ind}_fiber{self.fiber_ind}_amp{amp_ind}.dat'
        )
        with open(output_file_path, 'w') as activation_file:
            activation_file.write(f'{n_aps}')

    def save_variables(self, fiber: object, stimulation: ScaledStim, amp_ind: int = 0):
        """Write user-specified variables to file.

        :param fiber: instance of Fiber class
        :param stimulation: instance of ScaledStim class
        :param amp_ind: index of amplitude if protocol is FINITE_AMPLITUDES
        """
        # Put all recorded data into pandas DataFrame
        vm_data = pd.DataFrame([list(vm) for vm in fiber.vm if vm is not None])
        all_gating_data = [
            pd.DataFrame([list(g) for g in fiber.gating[gating_parameter] if g is not None])
            for gating_parameter in list(fiber.gating.keys())
        ]
        istim_data = pd.DataFrame(list(stimulation.istim_record))

        self.save_space_vm(amp_ind, fiber, vm_data)
        self.save_space_gating(all_gating_data, amp_ind, fiber)
        self.save_time_vm(amp_ind, stimulation, vm_data)
        self.save_time_gating(all_gating_data, amp_ind, stimulation)
        self.save_istim(amp_ind, istim_data, stimulation)
        self.save_aploctime(amp_ind, fiber)
        self.save_apendtimes(amp_ind, fiber)

    def space_header(self, header: list, var_type: str, dt: float, units: str = None):
        """Create header for F(x) files.

        :param header: list of string values to put at the header
        :param var_type: type of variable to be saved (Vm, h, mp, m, s)
        :param dt: float value of stimulation time step (ms)
        :param units: units of variable (mV, nA)
        """
        header.append('Node#')
        suffix = f'({units})' if units else ''
        for time in self.time_inds:
            header.append(f'{var_type}_time{int(time * dt)}ms{suffix}')

    def time_header(self, header: list, var_type: str, units: str = None):
        """Create header for F(t) files.

        :param header: list of string values to put at the header
        :param var_type: type of variable to be saved (Vm, h, mp, m, s)
        :param units: units of variable (mV, nA)
        """
        header.append('Time(ms)')
        suffix = f'({units})' if units else ''
        if var_type != 'istim':
            for node in self.node_inds:
                header.append(f'{var_type}_node{node + 1}{suffix}')
        elif var_type == 'istim':
            header.append(f'Istim({units})')

    def handle_header(self, save_type: str, var_type: str, dt: float = None, units: str = None):
        """Create a header for text file.

        :param save_type: function of variable to be saved, can be function of time ('time') or space ('space')
        :param var_type: type of variable to be saved (Vm, h, mp, m, s)
        :param units: units of variable (mV, nA)
        :param dt: float value for stimulation time step (ms)
        :return: list of column headers
        """
        header = []
        if save_type == 'space':  # F(x) - function of space
            self.space_header(header, var_type, dt, units)
        elif save_type == 'time':  # F(t) - function of time
            self.time_header(header, var_type, units)
        return header

    def save_istim(self, amp_ind: int, istim_data: list, stimulation: ScaledStim):
        """Save istim data to file.

        :param amp_ind: index of the amplitude in the finite amplitude protocol
        :param istim_data: list of istim data recorded during NEURON simulation list[nA]
        :param stimulation: instance of ScaledStim class
        """
        if self.istim:
            istim_path = os.path.join(
                self.output_path, f'Istim_inner{self.inner_ind}_fiber{self.fiber_ind}_amp{amp_ind}.dat'
            )
            istim_data.insert(0, 'Time', stimulation.time)
            istim_header = self.handle_header('time', 'istim', 'nA')
            istim_data.to_csv(istim_path, header=istim_header, sep='\t', float_format='%.6f', index=False)

    def save_time_gating(self, all_gating_data: list, amp_ind: int, stimulation: ScaledStim):
        """Save gating data as function of time to file.

        :param all_gating_data: list of lists containing float values for h, m, mp, and s gating parameters
        :param amp_ind: index of amplitude in finite amplitudes protocol
        :param stimulation: instance of ScaledStim class
        """
        if self.time_gating:
            gating_params = ['h', 'm', 'mp', 's']
            for gating_param, gating_data in zip(gating_params, all_gating_data):
                gating_time_path = os.path.join(
                    self.output_path,
                    f'gating_{gating_param}_time_inner{self.inner_ind}_fiber{self.fiber_ind}_amp{amp_ind}.dat',
                )
                gating_time_data = gating_data.T[self.node_inds]  # save data only at user-specified locations
                gating_time_data.insert(0, 'Time', stimulation.time)
                gating_time_header = self.handle_header('time', gating_param)
                gating_time_data.to_csv(
                    gating_time_path, header=gating_time_header, sep='\t', float_format='%.6f', index=False
                )

    def save_time_vm(self, amp_ind: int, stimulation: ScaledStim, vm_data: list):
        """Save membrane potential as a function of time to file.

        :param amp_ind: index of amplitude in finite amplitudes protocol
        :param stimulation: instance of ScaledStim class
        :param vm_data: list of float values for membrane potential as a function of time (mV)
        """
        if self.time_vm:
            vm_time_path = os.path.join(
                self.output_path, f'Vm_time_inner{self.inner_ind}_fiber{self.fiber_ind}_amp{amp_ind}.dat'
            )
            vm_time_data = vm_data.T[self.node_inds]  # save data only at user-specified locations
            vm_time_data.insert(0, 'Time', stimulation.time)
            vm_time_header = self.handle_header('time', 'vm', 'mV')
            vm_time_data.to_csv(vm_time_path, header=vm_time_header, sep='\t', float_format='%.6f', index=False)

    def save_space_gating(self, all_gating_data: list, amp_ind: int, fiber: _Fiber):
        """Save gating data as function of space to file.

        :param all_gating_data: list of lists containing float values for h, m, mp, and s gating parameters
        :param amp_ind: index of amplitude in finite amplitudes protocol
        :param fiber: instance of _Fiber class
        """
        if self.space_gating:
            gating_params = ['h', 'm', 'mp', 's']
            for gating_param, gating_data in zip(gating_params, all_gating_data):
                gating_space_path = os.path.join(
                    self.output_path,
                    f'gating_{gating_param}_space_inner{self.inner_ind}_fiber{self.fiber_ind}_amp{amp_ind}.dat',
                )
                gating_space_data = gating_data[self.time_inds]  # save data only at user-specified times
                if fiber.passive_end_nodes:
                    gating_space_data.insert(0, 'Node#', [*range(2, len(fiber.nodes))])
                else:
                    gating_space_data.insert(0, 'Node#', [*range(1, len(fiber.nodes) + 1)])
                gating_space_header = self.handle_header('space', gating_param)
                gating_space_data.to_csv(
                    gating_space_path, header=gating_space_header, sep='\t', float_format='%.6f', index=False
                )

    def save_space_vm(self, amp_ind, fiber, vm_data):
        """Save membrane potential data as function of space to file.

        :param amp_ind: index of amplitude in finite amplitudes protocol
        :param fiber: instance of _Fiber class
        :param vm_data: list of float values for membrane potential as a function of time (mV)
        """
        if self.space_vm:
            vm_space_path = os.path.join(
                self.output_path, f'Vm_space_inner{self.inner_ind}_fiber{self.fiber_ind}_amp{amp_ind}.dat'
            )
            vm_space_data = vm_data[self.time_inds]  # save data only at user-specified times
            if fiber.passive_end_nodes:
                vm_space_data.insert(0, 'Node#', [*range(2, len(fiber.nodes))])
            else:
                vm_space_data.insert(0, 'Node#', [*range(1, len(fiber.nodes) + 1)])
            vm_space_header = self.handle_header('space', 'vm', 'mV')
            vm_space_data.to_csv(vm_space_path, header=vm_space_header, sep='\t', float_format='%.6f', index=False)

    def save_aploctime(self, amp_ind: int, fiber: _Fiber):
        """Save time when last NoR AP was detected for each node to file.

        :param amp_ind: index of amplitude in finite amplitudes protocol
        :param fiber: instance of _Fiber class
        """
        if self.ap_loctime:
            aploctime_path = os.path.join(
                self.output_path, f'ap_loctime_inner{self.inner_ind}_fiber{self.fiber_ind}_amp{amp_ind}.dat'
            )
            with open(aploctime_path, 'w') as file:
                for node_ind, _ in enumerate(fiber.nodes):
                    if fiber.myelinated:
                        file.write(f"{fiber.apc[node_ind].time}")
                    else:
                        file.write(f"{fiber.apc[node_ind].time}")

    def save_apendtimes(self, amp_ind: int, fiber: _Fiber):
        """Save time that AP last propagated at two specified indices to file.

        :param amp_ind: index of amplitude in finite amplitudes protocol
        :param fiber: instance of _Fiber class
        """
        if self.ap_end_times:
            ap_end_times_path = os.path.join(
                self.output_path, f'Aptimes_inner{self.inner_ind}_fiber{self.fiber_ind}_amp{amp_ind}.dat'
            )
            with open(ap_end_times_path, 'w') as file:
                #
                file.write(f"{fiber.apc[self.ap_end_inds[0]].time} \t {fiber.apc[self.ap_end_inds[1]].time}")
