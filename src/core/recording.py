import os
import sys
from neuron import h

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '../')))
from src.utils import (Config, Configurable)
h.load_file('stdrun.hoc')

class Recording(Configurable):
    def __init__(self, fiber):
        self.time = h.Vector().record(h._ref_t)
        self.space = [i for i in range(1, fiber.axonnodes + 1)]
        self.vm = []

        self.gating_inds = [i for i in range(0, fiber.axonnodes)]
        if fiber.passive_end_nodes:
            del self.gating_inds[0]
            del self.gating_inds[-1]

        self.gating_h = [[] for i in self.gating_inds]
        self.gating_m = [[] for i in self.gating_inds]
        self.gating_mp = [[] for i in self.gating_inds]
        self.gating_s = [[] for i in self.gating_inds]
        self.gating = [self.gating_h, self.gating_m, self.gating_mp, self.gating_s]

        self.istim = []

        self.apc = []

    def record_ap(self, fiber):
        if fiber.fiber_type == 2:
            for i, node in enumerate(fiber.node):
                self.apc.append(h.APCount(node(0.5)))
                self.apc[i].thresh = fiber.search(Config.SIM, "protocol", "threshold", "value")
        else:
            for i, node in enumerate(fiber.sec):
                self.apc.append(h.APCount(node(0.5)))
                self.apc[i].thresh = fiber.search(Config.SIM, "protocol", "threshold", "value")


    def record_vm(self, fiber):
        for node_ind in range(0, fiber.axonnodes):
            if fiber.myelination:
                v_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_v)
                self.vm.append(v_node)
            else:
                v_node = h.Vector().record(fiber.sec[node_ind](0.5)._ref_v)
                self.vm.append(v_node)
        return

    def record_istim(self, i):
        self.istim.append(i)

    def record_gating(self, fiber, fix_passive=None):
        for j, node_ind in enumerate(self.gating_inds):
            self.gating_h[j].append(fiber.node[node_ind].hinf_axnodeMyel)
            self.gating_m[j].append(fiber.node[node_ind].minf_axnodeMyel)
            self.gating_mp[j].append(fiber.node[node_ind].mpinf_axnodeMyel)
            self.gating_s[j].append(fiber.node[node_ind].sinf_axnodeMyel)

        if fix_passive:
            for gating_vector in self.gating:
                n_vals = len(gating_vector[0])
                passive_gating_vals = [0 for i in range(0, n_vals)]
                gating_vector.insert(0, passive_gating_vals)
                gating_vector.append(passive_gating_vals)



    def ap_checker(self, fiber, find_block_thresh=False):
        ap_detect_location = fiber.search(Config.SIM, 'protocol', 'threshold', 'ap_detect_location')
        n_min_aps = fiber.search(Config.SIM, 'protocol', 'threshold', 'n_min_aps')
        node_index = int((fiber.axonnodes - 1) * ap_detect_location)

        if find_block_thresh:
            IntraStim_PulseTrain_delay = fiber.search(Config.SIM, 'intracellular_stim', 'times',
                                                     'IntraStim_PulseTrain_delay')
            if self.apc[node_index].time > IntraStim_PulseTrain_delay:
                ap_check = False
            else:
                ap_check = True
        else:
            if self.apc[node_index].n >= n_min_aps:
                ap_check = True
            else:
                ap_check = False
        return ap_check, self.apc[node_index].n