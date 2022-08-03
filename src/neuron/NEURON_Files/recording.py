from neuron import h

from src.utils.enums import Config
from src.utils.configurable import Configurable

h.load_file('stdrun.hoc')

class Recording(Configurable):
    def __init__(self, fiber):
        self.time = h.Vector().record(h._ref_t)
        self.space = [i for i in range(0, fiber.axonnodes)]
        self.vm = []

        self.gating_inds = [i for i in range(0, fiber.axonnodes)]
        if fiber.passive_end_nodes:
            del self.gating_inds[0]
            del self.gating_inds[-1]

        self.gating_h = []
        self.gating_m = []
        self.gating_mp = []
        self.gating_s = []
        self.gating = [self.gating_h, self.gating_m, self.gating_mp, self.gating_s]

        self.istim = []

        self.apc = []
        self.ap_end_count = []
        self.ap_end_times = []

    def reset(self):
        self.vm = []

        self.gating_h = []
        self.gating_m = []
        self.gating_mp = []
        self.gating_s = []
        self.gating = [self.gating_h, self.gating_m, self.gating_mp, self.gating_s]

        self.istim = []

        self.apc = []
        self.ap_end_count = []
        self.ap_end_times = []

    def record_ap(self, fiber):
        if fiber.myelination:
            for i, node in enumerate(fiber.node):
                self.apc.append(h.APCount(node(0.5)))
                thresh = fiber.search(Config.SIM, "protocol", "threshold", "value", optional=True)
                if thresh is not None:
                    self.apc[i].thresh = thresh
                else:
                    self.apc[i].thresh = -30
        else:
            for i, node in enumerate(fiber.sec):
                self.apc.append(h.APCount(node(0.5)))
                thresh = fiber.search(Config.SIM, "protocol", "threshold", "value", optional=True)
                if thresh is not None:
                    self.apc[i].thresh = thresh
                else:
                    self.apc[i].thresh = -30

    def record_ap_end_times(self, fiber, ap_end_inds, ap_end_thresh):
        self.ap_end_times = [h.Vector(), h.Vector()]
        for ap_end_vector, ap_end_ind in zip(self.ap_end_times, ap_end_inds):
            if fiber.myelination:
                ap_count = h.APCount(fiber.node[ap_end_ind](0.5))
                ap_count.thresh = ap_end_thresh
                ap_count.record(ap_end_vector)
                self.ap_end_count.append(ap_count)
            else:
                ap_end_min = h.APCount(fiber.sec[ap_end_ind](0.5))
                ap_end_min.thresh = ap_end_thresh
                ap_end_min.record(ap_end_vector)
                self.ap_end_count.append(ap_count)

    def record_vm(self, fiber):
        for node_ind in range(0, fiber.axonnodes):
            if fiber.myelination:
                v_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_v)
                self.vm.append(v_node)
            else:
                v_node = h.Vector().record(fiber.sec[node_ind](0.5)._ref_v)
                self.vm.append(v_node)
        return

    def record_istim(self, istim):
        self.istim = h.Vector().record(istim._ref_i)

    def record_gating(self, fiber, fix_passive=False):
        if fix_passive is False:
            for j, node_ind in enumerate(self.gating_inds):
                h_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_h_inf_axnode_myel)
                m_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_m_inf_axnode_myel)
                mp_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_mp_inf_axnode_myel)
                s_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_s_inf_axnode_myel)
                self.gating_h.append(h_node)
                self.gating_m.append(m_node)
                self.gating_mp.append(mp_node)
                self.gating_s.append(s_node)

        elif fix_passive and fiber.passive_end_nodes:
            for gating_vectors in self.gating:
                size = gating_vectors[0].size()
                passive_node = h.Vector(size, 0)
                gating_vectors.insert(0, passive_node)
                gating_vectors.append(passive_node)

    def ap_checker(self, fiber, find_block_thresh=False):
        ap_detect_location = fiber.search(Config.SIM, 'protocol', 'threshold', 'ap_detect_location', optional=True)
        if ap_detect_location is None:
            ap_detect_location = 0.9
        node_index = int((fiber.axonnodes - 1) * ap_detect_location)
        if find_block_thresh:
            IntraStim_PulseTrain_delay = fiber.search(Config.SIM, 'intracellular_stim', 'times',
                                                      'IntraStim_PulseTrain_delay')
            if self.apc[node_index].time > IntraStim_PulseTrain_delay:
                n_aps = False
            else:
                n_aps = True
        else:
            n_aps = self.apc[node_index].n
        return n_aps