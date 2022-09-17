"""The copyrights of this software are owned by Duke University.

Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""
from neuron import h
from src.utils import Config, Configurable

h.load_file('stdrun.hoc')


class Recording(Configurable):
    """Manage recording parameters for NEURON simulations."""

    def __init__(self, fiber: object):
        """Initialize Recording class.

        :param fiber: instance of Fiber class
        """
        self.time = h.Vector().record(h._ref_t)
        self.space = list(range(0, fiber.axonnodes))
        self.vm = []

        self.gating_inds = list(range(0, fiber.axonnodes))
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
        """Reset recording attributes in order to be used for subsequent runs."""
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

    def record_ap(self, fiber: object):
        # TODO: consider merging with record_ap_end_times
        """Create a list of NEURON APCount objects at all nodes along the axon.

        :param fiber: instance of Fiber class
        """
        # if fiber is myelinated, create APCount for each node of Ranvier
        if fiber.myelination:
            for i, node in enumerate(fiber.node):
                self.apc.append(h.APCount(node(0.5)))
                thresh = fiber.search(Config.SIM, "protocol", "threshold", "value", optional=True)
                if thresh is not None:
                    self.apc[i].thresh = thresh
                else:
                    # if thresh is not specified in sim.json (only allowed for FINITE_AMPLITUDES), use default value
                    self.apc[i].thresh = -30
        # if fiber is not myelinated (c-fiber model), create APCount for each segment
        else:
            for i, node in enumerate(fiber.sec):
                self.apc.append(h.APCount(node(0.5)))
                thresh = fiber.search(Config.SIM, "protocol", "threshold", "value", optional=True)
                if thresh is not None:
                    self.apc[i].thresh = thresh
                else:
                    # if thresh is not specified in sim.json (only allowed for FINITE_AMPLITUDES), use default value
                    self.apc[i].thresh = -30

    def record_ap_end_times(self, fiber: object, ap_end_inds: list, ap_end_thresh: float):
        """Record when action potential occurs at specified indices. For 'end_ap_times' in sim.json.

        :param fiber: instance of Fiber class
        :param ap_end_inds: list of user-specified indices to record APs
        :param ap_end_thresh: threshold value for action potentials
        """
        # Create vectors to save ap times to
        self.ap_end_times = [h.Vector(), h.Vector()]
        for ap_end_vector, ap_end_ind in zip(self.ap_end_times, ap_end_inds):
            if fiber.myelination:
                # if myelinated, create APCount at node of Ranvier
                ap_count = h.APCount(fiber.node[ap_end_ind](0.5))
                ap_count.thresh = ap_end_thresh
                ap_count.record(ap_end_vector)  # save AP times detected by APCount to vector
                self.ap_end_count.append(ap_count)
            else:
                # if unmyelinated, create APCount at axon segment
                ap_end_min = h.APCount(fiber.sec[ap_end_ind](0.5))
                ap_end_min.thresh = ap_end_thresh
                ap_end_min.record(ap_end_vector)  # save AP times detected by APCount to vector
                self.ap_end_count.append(ap_count)

    def record_vm(self, fiber):
        """Record membrane voltage (mV) along the axon.

        :param fiber: instance of Fiber class
        """
        for node_ind in range(0, fiber.axonnodes):
            if fiber.myelination:
                v_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_v)
                self.vm.append(v_node)
            else:
                v_node = h.Vector().record(fiber.sec[node_ind](0.5)._ref_v)
                self.vm.append(v_node)
        return

    def record_istim(self, istim: object):
        """Record applied intracellular stimulation (nA).

        :param istim: instance of intracellular stimulation object
        """
        self.istim = h.Vector().record(istim._ref_i)

    def record_gating(self, fiber: object, fix_passive: bool = False):
        """Record gating parameters (h, m, mp, s) for myelinated fiber types.

        :param fiber: instance of Fiber class
        :param fix_passive: true if fiber has passive end nodes, false otherwise
        """
        if fix_passive is False:
            # Set up recording vectors for h, m, mp, and s gating parameters all along the axon
            for node_ind in self.gating_inds:
                h_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_h_inf_axnode_myel)
                m_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_m_inf_axnode_myel)
                mp_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_mp_inf_axnode_myel)
                s_node = h.Vector().record(fiber.node[node_ind](0.5)._ref_s_inf_axnode_myel)
                self.gating_h.append(h_node)
                self.gating_m.append(m_node)
                self.gating_mp.append(mp_node)
                self.gating_s.append(s_node)

        # If fiber has passive end nodes, then insert vectors of 0's to output
        elif fix_passive and fiber.passive_end_nodes:
            for gating_vectors in self.gating:
                size = gating_vectors[0].size()
                passive_node = h.Vector(size, 0)
                gating_vectors.insert(0, passive_node)
                gating_vectors.append(passive_node)

    def ap_checker(self, fiber: object, find_block_thresh: bool = False) -> int:
        """Check to see if an action potential occurred at the end of a run.

        :param fiber: instance of Fiber class
        :param find_block_thresh: true if BLOCK_THRESHOLD protocol, false otherwise
        :return: number of action potentials that occurred
        """
        # Determine user-specified location along axon to check for action potential
        ap_detect_location = fiber.search(Config.SIM, 'protocol', 'threshold', 'ap_detect_location', optional=True)
        if ap_detect_location is None:
            ap_detect_location = 0.9
        node_index = int((fiber.axonnodes - 1) * ap_detect_location)

        if find_block_thresh:
            intrastim_pulsetrain_delay = fiber.search(
                Config.SIM, 'intracellular_stim', 'times', 'IntraStim_PulseTrain_delay'
            )
            if self.apc[node_index].time > intrastim_pulsetrain_delay:
                n_aps = 0  # False - block did not occur
            else:
                n_aps = 1  # True - block did occur
        else:
            n_aps = self.apc[node_index].n
        return n_aps
