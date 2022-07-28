import os
import sys
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '../../')))
from src.utils import (Config, Configurable)

class Saving(Configurable):
    def __init__(self):
        Configurable.__init__(self)

        self.space_vm = None
        self.space_gating = None
        self.time_inds = []
        self.time_vm = None
        self.time_gating = None
        self.istim = None
        self.locs = []
        self.node_inds = []
        self.ap_end_times = None
        self.ap_end_inds = []
        self.ap_end_thresh = None
        self.runtime = None
        self.output_path = None
        return

    def inherit(self, sim_path, dt, fiber):
        self.space_vm = fiber.search(Config.SIM, "saving", "space", "vm")
        self.space_gating = fiber.search(Config.SIM, "saving", "space", "gating")
        sim_times = fiber.search(Config.SIM, "saving", "space", "times")
        self.time_inds = [int(t/dt) for t in sim_times]
        self.time_inds.sort()
        self.time_vm = fiber.search(Config.SIM, "saving", "time", "vm")
        self.time_gating = fiber.search(Config.SIM, "saving", "time", "gating")
        self.istim = fiber.search(Config.SIM, "saving", "time", "istim")
        self.locs = fiber.search(Config.SIM, "saving", "time", "locs")
        if self.locs != 'all':
            self.node_inds = [int((fiber.axonnodes-1)*loc) for loc in self.locs]
        elif self.locs == 'all':
            self.node_inds = [i for i in range(0, fiber.axonnodes)]
        self.node_inds.sort()
        saving = fiber.search(Config.SIM, "saving")
        if saving.get("end_ap_times") is not None:
            self.ap_end_times = True
            loc_min = fiber.search(Config.SIM, "saving", "end_ap_times", "loc_min")
            loc_max = fiber.search(Config.SIM, "saving", "end_ap_times", "loc_max")
            node_ind_min = int((fiber.axonnodes-1)*loc_min)
            node_ind_max = int((fiber.axonnodes - 1)*loc_max)
            self.ap_end_inds = [node_ind_min, node_ind_max]
            self.ap_end_thresh = fiber.search(Config.SIM, "saving", "end_ap_times", "threshold")
        self.runtime = fiber.search(Config.SIM, "saving", "runtimes")
        self.output_path = os.path.join(sim_path, 'data', 'outputs')
        return

    def saveThresh(self, fiber, thresh):
        # save threshold to submit/n_sims/#/data/outputs/thresh_inner#_fiber#.dat
        thresh_path = os.path.join(self.output_path, 'thresh_inner{}_fiber{}.dat'.format(fiber.inner_ind, fiber.fiber_ind))
        thresh_file = open(thresh_path, 'w')
        thresh_file.write("{:.6f} mA".format(thresh))
        thresh_file.close()

    def saveRuntime(self, fiber, runtime, amp_ind=0):
        if self.runtime:
            runtimes_path = os.path.join(self.output_path, 'runtime_inner{}_fiber{}_amp{}.dat'.format(fiber.inner_ind, fiber.fiber_ind, amp_ind))
            file = open(runtimes_path, 'w')
            file.write('{:.3f}s'.format(runtime))
            file.close()


    def write2file(self, recording, fiber, dt, amp_ind=0):
        def create_header(save_type: str, var_type: str, units: str = None):
            header = []
            if save_type == 'space':  # F(x) - function of space
                header.append('Node#')
                for time in self.time_inds:
                    if units is not None:
                        header.append(f'{var_type}_time{int(time * dt)}ms({units})')
                    else:
                        header.append(f'{var_type}_time{int(time * dt)}ms')
            elif save_type == 'time':  # F(t) - function of time
                header.append('Time(ms)')
                if var_type != 'istim':
                    for node in self.node_inds:
                        if units is not None:
                            header.append(f'{var_type}_node{node + 1}({units})')
                        else:
                            header.append(f'{var_type}_node{node + 1}')
                elif var_type == 'istim':
                    header.append(f'Istim({units})')
            return header

        vm_data = pd.DataFrame(recording.vm)
        all_gating_data = [pd.DataFrame(gating_vector) for gating_vector in recording.gating]
        istim_data = pd.DataFrame(recording.istim)

        if self.space_vm:
            vm_space_path = os.path.join(self.output_path, 'Vm_space_inner{}_fiber{}_amp{}.dat'.format(fiber.inner_ind, fiber.fiber_ind, amp_ind))
            vm_space_data = vm_data[self.time_inds]
            vm_space_data.insert(0, 'Node#', recording.space)
            vm_space_header = create_header('space', 'vm', 'mV')
            vm_space_data.to_csv(vm_space_path, header=vm_space_header, sep='\t', float_format='%.6f', index=False)

        if self.space_gating:
            gating_params = ['h', 'm', 'mp', 's']
            for gating_param, gating_data in zip(gating_params, all_gating_data):
                gating_space_path = os.path.join(self.output_path, 'gating_{}_space_inner{}_fiber{}_amp{}.dat'.format(gating_param, fiber.inner_ind, fiber.fiber_ind, amp_ind))
                gating_space_data = gating_data[self.time_inds]
                gating_space_data.insert(0, 'Node#', recording.space)
                gating_space_header = create_header('space', gating_param)
                gating_space_data.to_csv(gating_space_path, header=gating_space_header, sep='\t', float_format='%.6f', index=False)

        if self.time_vm:
            vm_time_path = os.path.join(self.output_path, 'Vm_time_inner{}_fiber{}_amp{}.dat'.format(fiber.inner_ind, fiber.fiber_ind, amp_ind))
            vm_time_data = vm_data.T[self.node_inds]
            vm_time_data.insert(0, 'Time', recording.time)
            vm_time_header = create_header('time', 'vm', 'mV')
            vm_time_data.to_csv(vm_time_path, header=vm_time_header, sep='\t', float_format='%.6f', index=False)

        if self.time_gating:
            gating_params = ['h', 'm', 'mp', 's']
            for gating_param, gating_data in zip(gating_params, all_gating_data):
                gating_time_path = os.path.join(self.output_path,
                                                 'gating_{}_time_inner{}_fiber{}_amp{}.dat'.format(gating_param, fiber.inner_ind,
                                                                                              fiber.fiber_ind, amp_ind))
                gating_time_data = gating_data.T[self.node_inds]
                gating_time_data.insert(0, 'Time', recording.time)
                gating_time_header = create_header('time', gating_param)
                gating_time_data.to_csv(gating_time_path, header=gating_time_header, sep='\t', float_format='%.6f', index=False)

        if self.istim:
            istim_path = os.path.join(self.output_path, 'Istim_inner{}_fiber{}_amp{}.dat'.format(fiber.inner_ind, fiber.fiber_ind, amp_ind))
            istim_data.insert(0, 'Time', recording.time)
            istim_header = create_header('time', 'istim', 'nA')
            istim_data.to_csv(istim_path, header=istim_header, sep='\t', float_format='%.6f', index=False)

        if self.ap_end_times:
            ap_end_times_path = os.path.join(self.output_path, 'Aptimes_inner{}_fiber{}_amp{}.dat'.format(fiber.inner_ind, fiber.fiber_ind, amp_ind))
            file = open(ap_end_times_path, 'w')
            file.write(f'Node{self.ap_end_inds[0] + 1} \t Node{self.ap_end_inds[1] + 1} \n')

            min_size, max_size = recording.ap_end_times[0].size(), recording.ap_end_times[1].size()
            n_rows = max(min_size, max_size)
            for i in range(0, n_rows):
                if i < min_size and i < max_size:
                    file.write('{:.6f} \t {:.6f} \n'.format(recording.ap_end_times[0][i], recording.ap_end_times[1][i]))
                elif min_size > i >= max_size:
                    file.write('{:.6f} \t  nan \n'.format(recording.ap_end_times[0][i]))
                elif min_size <= i < max_size:
                    file.write('nan \t {:.6f} \n'.format(recording.ap_end_times[1][i]))
            file.close()