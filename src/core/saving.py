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
        self.time_Istim = None
        self.locs = []
        self.node_inds = []
        self.runtime = None
        self.time_vector = []
        self.vm_vectors = []

    def inherit(self, sim_path, fiber_axonnodes, fiber_dz):
        self.space_vm = self.search(Config.SIM, "saving", "space", "vm")
        self.space_gating = self.search(Config.SIM, "saving", "space", "gating")
        self.times = self.search(Config.SIM, "saving", "space", "times")
        self.time_vm = self.search(Config.SIM, "saving", "time", "vm")
        self.time_gating = self.search(Config.SIM, "saving", "time", "gating")
        self.time_Istim = self.search(Config.SIM, "saving", "time", "istim")
        self.locs = self.search(Config.SIM, "saving", "time", "locs")
        for loc in self.locs:
            self.node_inds.append(int((fiber_axonnodes-1)*loc))
        self.runtime = self.search(Config.SIM, "saving", "runtimes")
        self.output_path = os.path.join(sim_path, 'data', 'outputs')

    def save_data(self, save_type: str, var_type: str, x_vector: list, y_vectors: list):
        if save_type == 'space':
            pass
        elif save_type == 'time':
            self.time_vector = list(x_vector)
            if var_type == 'vm':
                self.vm_vectors = [list(x) for x in y_vectors]
        return

    def write2file(self, inner_ind, fiber_ind, amp_ind=0):
        def create_header(save_type: str, var_type: str, units: str):
            if save_type == 'space':  # F(x) - function of space
                header = 'Node# \t '
                for time in self.times:
                    header = header + '{}_time{}ms ({}) \t '.format(var_type, time, units)

            elif save_type == 'time':  # F(t) - function of time
                header = 'Time(ms) \t '
                for node in self.node_inds:  # todo: turn locs into self.index
                    header = header + '{}_node{} ({}) \t '.format(var_type, node, units)

            return header + '\n'

        def write_lines(file_path, lines):
            file = open(file_path, 'w')
            for line in lines:
                file.write(line)
            file.close()

        if self.space_vm:
            vm_space_path = os.path.join(self.output_path, 'Vm_space_inner{}_fiber{}_amp{}.dat'.format(inner_ind, fiber_ind, amp_ind))
        if self.space_gating:
            gating_space_path = os.path.join(self.output_path, 'gating_{}_space_inner{}_fiber{}_amp{}.dat'.format(gating_param, inner_ind, fiber_ind, amp_ind))
        if self.time_vm:
            vm_time_path = os.path.join(self.output_path, 'Vm_time_inner{}_fiber{}_amp{}.dat'.format(inner_ind, fiber_ind, amp_ind))
            lines = []  # list of strings with lines to write to file
            lines.append(create_header('time', 'vm', 'nA'))
            for t, time in enumerate(self.time_vector):
                line = '{:.6f} \t '.format(time)
                for vm_node in self.vm_vectors:
                    line = line + '{:.6f} \t '.format(vm_node[t])
                lines.append(line + '\n')
            write_lines(vm_time_path, lines)

        if self.time_gating:
            gating_time_path = os.path.join(self.output_path,
                                             'gating_{}_time_inner{}_fiber{}_amp{}.dat'.format(gating_param, inner_ind,
                                                                                          fiber_ind, amp_ind))
        if self.time_Istim:
            istim_path = os.path.join(self.output_path, 'Istim_inner{}_fiber{}_amp{}.dat'.format(inner_ind, fiber_ind, amp_ind))
        if self.runtime:
            runtimes_path = os.path.join(self.output_path, 'runtime_inner{}_fiber{}_amp{}.dat'.format(inner_ind, fiber_ind, amp_ind))
