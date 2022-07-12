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
        self.time_vector = []
        self.time_vm_vectors = []
        self.time_gating_h_vectors =[]
        self.time_gating_m_vectors = []
        self.time_gating_mp_vectors = []
        self.time_gating_s_vectors = []
        self.space_gating_h_vectors =[]
        self.space_gating_m_vectors = []
        self.space_gating_mp_vectors = []
        self.space_gating_s_vectors = []
        self.space_vector = []
        self.space_vm_vectors = []

    def inherit(self, sim_path, dt, fiber_axonnodes):
        self.space_vm = self.search(Config.SIM, "saving", "space", "vm")
        self.space_gating = self.search(Config.SIM, "saving", "space", "gating")
        self.times = self.search(Config.SIM, "saving", "space", "times")
        self.times = [int(t/dt) for t in self.times]
        self.times.sort()
        self.time_vm = self.search(Config.SIM, "saving", "time", "vm")
        self.time_gating = self.search(Config.SIM, "saving", "time", "gating")
        self.istim = self.search(Config.SIM, "saving", "time", "istim")
        self.locs = self.search(Config.SIM, "saving", "time", "locs")
        if self.locs != 'all':
            self.node_inds = [int((fiber_axonnodes-1)*loc) for loc in self.locs]
        elif self.locs == 'all':
            self.node_inds = [i for i in range(0, fiber_axonnodes)]
        self.node_inds.sort()
        self.runtime = self.search(Config.SIM, "saving", "runtimes")
        self.output_path = os.path.join(sim_path, 'data', 'outputs')

    def save_data(self, save_type: str, var_type: str, x_vector: list, y_vectors: list):
        if save_type == 'space':
            if len(self.space_vector) == 0:
                self.space_vector = x_vector
            if var_type == 'vm':
                self.space_vm_vectors = y_vectors
            if var_type == 'gating_h':
                self.space_gating_h_vectors = y_vectors
            if var_type == 'gating_m':
                self.space_gating_m_vectors = y_vectors
            if var_type == 'gating_mp':
                self.space_gating_mp_vectors = y_vectors
            if var_type == 'gating_s':
                self.space_gating_s_vectors = y_vectors
        elif save_type == 'time':
            if len(self.time_vector) == 0:
                self.time_vector = x_vector
            if var_type == 'vm':
                self.time_vm_vectors = y_vectors
            if var_type == 'gating_h':
                self.time_gating_h_vectors = y_vectors
            if var_type == 'gating_m':
                self.time_gating_m_vectors = y_vectors
            if var_type == 'gating_mp':
                self.time_gating_mp_vectors = y_vectors
            if var_type == 'gating_s':
                self.time_gating_s_vectors = y_vectors
        return

    def write2file(self, runtime, inner_ind, fiber_ind, dt, amp_ind=0):
        def create_header(save_type: str, var_type: str, units: str = None):
            if save_type == 'space':  # F(x) - function of space
                header = 'Node# \t '
                for time in self.times:
                    if units is not None:
                        header = header + '{}_time{}ms({}) \t '.format(var_type, time*dt, units)
                    else:
                        header = header + '{}_time{}ms \t '.format(var_type, time*dt)

            elif save_type == 'time':  # F(t) - function of time
                header = 'Time(ms) \t '
                for node in self.node_inds:
                    if units is not None:
                        header = header + '{}_node{}({}) \t '.format(var_type, node + 1, units)
                    else:
                        header = header + '{}_node{} \t '.format(var_type, node + 1)
            return header + '\n'

        def write_lines(file_path, lines):
            file = open(file_path, 'w')
            for line in lines:
                file.write(line)
            file.close()

        if self.space_vm:
            vm_space_path = os.path.join(self.output_path, 'Vm_space_inner{}_fiber{}_amp{}.dat'.format(inner_ind, fiber_ind, amp_ind))
            lines = []  # list of strings with lines to write to file
            lines.append(create_header('space', 'vm', units='mV'))
            for t, space in enumerate(self.space_vector):
                line = '{} \t \t '.format(space)
                for vm_time in self.space_vm_vectors:
                    line = line + '{:.6f} \t '.format(vm_time[t])
                lines.append(line + '\n')
            write_lines(vm_space_path, lines)

        if self.space_gating:
            gating_params, gating_vectors = ['h', 'm', 'mp', 's'], [self.space_gating_h_vectors,
                                                                    self.space_gating_m_vectors,
                                                                    self.space_gating_mp_vectors,
                                                                    self.space_gating_s_vectors]
            for gating_param, gating_vector in zip(gating_params, gating_vectors):
                gating_space_path = os.path.join(self.output_path, 'gating_{}_space_inner{}_fiber{}_amp{}.dat'.format(gating_param, inner_ind, fiber_ind, amp_ind))
                lines = []  # list of strings with lines to write to file
                lines.append(create_header('space', gating_param))
                for t, space in enumerate(self.space_vector):
                    line = '{} \t \t '.format(space)
                    for gating_time in gating_vector:
                        line = line + '{:.6f} \t '.format(gating_time[t])
                    lines.append(line + '\n')
                write_lines(gating_space_path, lines)
        if self.time_vm:
            vm_time_path = os.path.join(self.output_path, 'Vm_time_inner{}_fiber{}_amp{}.dat'.format(inner_ind, fiber_ind, amp_ind))
            lines = []  # list of strings with lines to write to file
            lines.append(create_header('time', 'vm', units='mV'))
            for t, time in enumerate(self.time_vector):
                line = '{:.6f} \t '.format(time)
                for vm_node in self.time_vm_vectors:
                    line = line + '{:.6f} \t '.format(vm_node[t])
                lines.append(line + '\n')
            write_lines(vm_time_path, lines)

        if self.time_gating:
            gating_params, gating_vectors = ['h', 'm', 'mp', 's'], [self.time_gating_h_vectors, self.time_gating_m_vectors,
                                                                    self.time_gating_mp_vectors, self.time_gating_s_vectors]
            for gating_param, gating_vector in zip(gating_params, gating_vectors):
                gating_time_path = os.path.join(self.output_path,
                                                 'gating_{}_time_inner{}_fiber{}_amp{}.dat'.format(gating_param, inner_ind,
                                                                                              fiber_ind, amp_ind))
                gating_vector = [vec for i, vec in enumerate(gating_vector) if i in self.node_inds]
                lines = []  # list of strings with lines to write to file
                lines.append(create_header('time', gating_param))
                for t, time in enumerate(self.time_vector):
                    line = '{:.6f} \t '.format(time)
                    for gating_node in gating_vector:
                        line = line + '{:.6f} \t '.format(gating_node[t])
                    lines.append(line + '\n')
                write_lines(gating_time_path, lines)

        if self.istim:
            istim_path = os.path.join(self.output_path, 'Istim_inner{}_fiber{}_amp{}.dat'.format(inner_ind, fiber_ind, amp_ind))
        if self.runtime:
            runtimes_path = os.path.join(self.output_path, 'runtime_inner{}_fiber{}_amp{}.dat'.format(inner_ind, fiber_ind, amp_ind))
            file = open(runtimes_path, 'w')
            file.write('{:.3f}s'.format(runtime))
            file.close()
