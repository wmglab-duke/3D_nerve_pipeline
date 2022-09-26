"""Defines Stimulation class.

The copyrights of this software are owned by Duke University. Please
refer to the LICENSE and README.md files for licensing instructions. The
source code can be found on the following GitHub repository:
https://github.com/wmglab-duke/ascent
"""

from neuron import h
from src.utils import Config, Configurable

h.load_file('stdrun.hoc')


class Stimulation(Configurable):
    """Manage stimulation of NEURON simulations."""

    def __init__(self):
        """Initialize Stimulation class."""
        # initialize Configurable super class
        Configurable.__init__(self)

        self.potentials = []
        self.waveform = []
        self.dt = None
        self.tstop = None
        self.istim = None
        return

    def load_potentials(self, potentials_path: str):
        # TODO: have run_controls.py load the potentials and pass them in as a numpy array
        """Create Ve(x) -- vector of potentials from FEM.

        :param potentials_path: file name containing Extracellular Stim potentials data
        :return: Stimulation object
        :raises Exception: size of the potentials file does not match the Fiber class
        """
        with open(potentials_path, 'r') as potentials_file:
            axontotal = int(potentials_file.readline())
            file_lines = potentials_file.read().splitlines()
            self.potentials = [float(i) * 1000 for i in file_lines]  # Need to convert to V -> mV

        if len(self.potentials) != axontotal:
            raise Exception("Need axontotal from potentials file to match axontotal used in Python")
        return self

    def load_waveform(self, waveform_path: str):
        # TODO: have run_controls.py load the waveform and pass them in as a numpy array
        """Create I(t) -- vector of amplitudes at each time step of the FEM. Read in simulation time step and time stop.

        :param waveform_path: file name containing Extracellular Stim waveform data
        :return: Stimulation object
        """
        with open(waveform_path, 'r') as waveform_file:
            self.dt = float(waveform_file.readline().strip())  # time step
            self.tstop = int(waveform_file.readline().strip())  # stop time
            file_lines = waveform_file.read().splitlines()
            self.waveform = [float(i) for i in file_lines]
        return self

    def apply_intracellular(self, fiber: object):
        """Create instance of trainIClamp for intracellular stimulation.

        :param fiber: instance of Fiber class
        :return: instance of Stimulation class
        """
        intrastim_pulsetrain_ind = fiber.search(Config.SIM, 'intracellular_stim', 'ind')
        if fiber.myelination:
            ind = intrastim_pulsetrain_ind*11
        else:
            ind = intrastim_pulsetrain_ind
        intracellular_stim = h.trainIClamp(fiber.segments[ind](0.5))
        intracellular_stim.delay = fiber.search(Config.SIM, 'intracellular_stim', 'times', 'IntraStim_PulseTrain_delay')
        intracellular_stim.PW = fiber.search(Config.SIM, 'intracellular_stim', 'times', 'pw')
        intracellular_stim.train = fiber.search(Config.SIM, 'intracellular_stim', 'times', 'IntraStim_PulseTrain_dur')
        intracellular_stim.freq = fiber.search(Config.SIM, 'intracellular_stim', 'pulse_repetition_freq')
        intracellular_stim.amp = fiber.search(Config.SIM, 'intracellular_stim', 'amp')
        self.istim = intracellular_stim
        return self

    def initialize_extracellular(self, fiber: object):
        """Set extracellular stimulation values to zero along entire fiber.

        :param fiber: instance of Fiber class
        """
        for segment in fiber.segments:
            segment(0.5).e_extracellular = 0
        return

    def update_extracellular(self, fiber: object, e_stims: str):
        """Update the applied extracellular stimulation all along the fiber length.

        :param fiber: instance of Fiber class
        :param e_stims: list of extracellular stimulations to apply along fiber length
        """
        for x, segment in enumerate(fiber.segments):
            segment(0.5).e_extracellular = e_stims[x]
