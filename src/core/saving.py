import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '../')))
from src.utils import (Config, Configurable)

class Saving(Configurable):
    def __init__(self):
        Configurable.__init__(self)

        self.space_vm = None
        self.space_gating = None
        self.times = []
        self.time_vm = None
        self.time_gating = None
        self.istim = None
        self.locs = []
        self.node_inds = []
        self.runtime = None
        return

    def inherit(self, sim_path, dt, fiber):
        self.space_vm = fiber.search(Config.SIM, "saving", "space", "vm")
        self.space_gating = fiber.search(Config.SIM, "saving", "space", "gating")
        self.times = fiber.search(Config.SIM, "saving", "space", "times")
        self.times = [int(t/dt) for t in self.times]
        self.times.sort()
        self.time_vm = fiber.search(Config.SIM, "saving", "time", "vm")
        self.time_gating = fiber.search(Config.SIM, "saving", "time", "gating")
        self.istim = fiber.search(Config.SIM, "saving", "time", "istim")
        self.locs = fiber.search(Config.SIM, "saving", "time", "locs")
        if self.locs != 'all':
            self.node_inds = [int((fiber.axonnodes-1)*loc) for loc in self.locs]
        elif self.locs == 'all':
            self.node_inds = [i for i in range(0, fiber.axonnodes)]
        self.node_inds.sort()
        self.runtime = fiber.search(Config.SIM, "saving", "runtimes")
        self.output_path = os.path.join(sim_path, 'data', 'outputs')
        return

    def saveThresh(self, fiber, thresh):
        # save threshold to submit/n_sims/#/data/outputs/thresh_inner#_fiber#.dat
        thresh_path = os.path.join(self.output_path, 'thresh_inner{}_fiber{}.dat'.format(fiber.inner_ind, fiber.fiber_ind))
        thresh_file = open(thresh_path, 'w')
        thresh_file.write("{:.6f} mA".format(thresh))
        thresh_file.close()

    def write2file(self, recording, runtime, fiber, dt, amp_ind=0):
        def create_header(save_type: str, var_type: str, units: str = None):
            if save_type == 'space':  # F(x) - function of space
                header = 'Node# \t '
                for time in self.times:
                    if units is not None:
                        header = header + '{}_time{}ms({}) \t '.format(var_type, int(time*dt), units)
                    else:
                        header = header + '{}_time{}ms \t '.format(var_type, int(time*dt))
            elif save_type == 'time':  # F(t) - function of time
                header = 'Time(ms) \t '
                if var_type != 'istim':
                    for node in self.node_inds:
                        if units is not None:
                            header = header + '{}_node{}({}) \t '.format(var_type, node + 1, units)
                        else:
                            header = header + '{}_node{} \t '.format(var_type, node + 1)
                elif var_type == 'istim':
                    header = header + 'Istim({})'.format(units)
            return header + '\n'
        
        def create_file_lines(x_vector: list, y_vector: list, save_type: str, var_type: str, units: str = None):
            file_lines = []  # list of strings with lines to write to file
            file_lines.append(create_header(save_type, var_type, units))
            for t, x_val in enumerate(x_vector):
                if save_type == 'space': line = '{} \t \t '.format(x_val)
                elif save_type == 'time': line = '{:.6f} \t \t '.format(x_val)
                if var_type != 'istim':
                    for y_val in y_vector:
                        line = line + '{:.6f} \t '.format(y_val[t])
                else: line = line + '{:.6f} \t '.format(y_vector[t])
                file_lines.append(line + '\n')
            return file_lines

        def write_file(file_path, lines):
            file = open(file_path, 'w')
            for line in lines:
                file.write(line)
            file.close()

        def create_space_vectors(saving_times, recording_vectors):
            space_vectors = [[] for t in self.times]
            for space_vec, time in zip(space_vectors, saving_times):
                for vec in recording_vectors:
                    space_vec.append(vec[time])
            return space_vectors

        def create_time_vectors(saving_inds, recording_vectors):
            return [vec for i, vec in enumerate(recording_vectors) if i in saving_inds]

        if self.space_vm:
            vm_space_path = os.path.join(self.output_path, 'Vm_space_inner{}_fiber{}_amp{}.dat'.format(fiber.inner_ind, fiber.fiber_ind, amp_ind))
            vm_space_vectors = create_space_vectors(self.times, recording.vm)
            file_lines = create_file_lines(recording.space, vm_space_vectors, 'space', 'vm', 'mV')
            write_file(vm_space_path, file_lines)

        if self.space_gating:
            gating_params = ['h', 'm', 'mp', 's']
            for gating_param, gating_vector in zip(gating_params, recording.gating):
                gating_space_path = os.path.join(self.output_path, 'gating_{}_space_inner{}_fiber{}_amp{}.dat'.format(gating_param, fiber.inner_ind, fiber.fiber_ind, amp_ind))
                gating_space_vectors = create_space_vectors(self.times, gating_vector)
                file_lines = create_file_lines(recording.space, gating_space_vectors, 'space', gating_param)
                write_file(gating_space_path, file_lines)

        if self.time_vm:
            vm_time_path = os.path.join(self.output_path, 'Vm_time_inner{}_fiber{}_amp{}.dat'.format(fiber.inner_ind, fiber.fiber_ind, amp_ind))
            vm_time_vectors = create_time_vectors(self.node_inds, recording.vm)
            time = list(recording.time)
            file_lines = create_file_lines(time, vm_time_vectors, 'time', 'vm', 'mV')
            write_file(vm_time_path, file_lines)

        if self.time_gating:
            gating_params = ['h', 'm', 'mp', 's']
            for gating_param, gating_vector in zip(gating_params, recording.gating):
                gating_time_path = os.path.join(self.output_path,
                                                 'gating_{}_time_inner{}_fiber{}_amp{}.dat'.format(gating_param, fiber.inner_ind,
                                                                                              fiber.fiber_ind, amp_ind))
                gating_time_vectors = create_time_vectors(self.node_inds, gating_vector)
                file_lines = create_file_lines(time, gating_time_vectors, 'time', gating_param)
                write_file(gating_time_path, file_lines)

        if self.istim:
            istim_path = os.path.join(self.output_path, 'Istim_inner{}_fiber{}_amp{}.dat'.format(fiber.inner_ind, fiber.fiber_ind, amp_ind))
            file_lines = create_file_lines(time, recording.istim, 'time', 'istim', 'nA')
            write_file(istim_path, file_lines)

        if self.runtime:
            runtimes_path = os.path.join(self.output_path, 'runtime_inner{}_fiber{}_amp{}.dat'.format(fiber.inner_ind, fiber.fiber_ind, amp_ind))
            file = open(runtimes_path, 'w')
            file.write('{:.3f}s'.format(runtime))
            file.close()
