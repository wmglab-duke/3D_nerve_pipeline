from neuron import h
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '../')))
from src.utils import (Config, Configurable)
h.load_file('stdrun.hoc')

class Stimulation(Configurable):
    def __init__(self):
        Configurable.__init__(self)

        self.potentials = []
        self.waveform = []
        self.dt = None
        self.tstop = None
        self.istim = None
        return

    def load_potentials(self, potentials_path):
        """
        :param potentials_path: file name containing Extracellular Stim potentials data
        Creates Ve(x) -- vector of potentials from FEM
        """
        potentials_file = open(potentials_path, 'r')
        axontotal = int(potentials_file.readline())
        file_lines = potentials_file.read().splitlines()
        self.potentials = [float(i)*1000 for i in file_lines]  # Need to convert to V -> mV
        potentials_file.close()

        if len(self.potentials) != axontotal:
            raise Exception("Need axontotal from potentials file to match axontotal used in Python")
        return self

    def load_waveform(self, waveform_path):
        """
        :param waveform_path: file name containing Extracellular Stim waveform data
        Creates I(t) -- vector of amplitudes at each time step of the FEM
        Also reads in time step and time stop
        """
        waveform_file = open(waveform_path, 'r')
        self.dt = float(waveform_file.readline().strip())           # time step
        self.tstop = int(waveform_file.readline().strip())          # stop time
        file_lines = waveform_file.read().splitlines()
        self.waveform = [float(i) for i in file_lines]
        waveform_file.close()
        return self

    def apply_intracellular(self, fiber):
        if fiber.myelination:
            fiber_sections = fiber.node
        else:
            fiber_sections = fiber.sec
        IntraStim_PulseTrain_ind = fiber.search(Config.SIM, 'intracellular_stim', 'ind')
        intracellular_stim = h.trainIClamp(fiber_sections[IntraStim_PulseTrain_ind](0.5))
        intracellular_stim.delay = fiber.search(Config.SIM, 'intracellular_stim', 'times', 'IntraStim_PulseTrain_delay')
        intracellular_stim.PW = fiber.search(Config.SIM, 'intracellular_stim', 'times', 'pw')
        intracellular_stim.train = fiber.search(Config.SIM, 'intracellular_stim', 'times', 'IntraStim_PulseTrain_dur')
        intracellular_stim.freq = fiber.search(Config.SIM, 'intracellular_stim', 'pulse_repetition_freq')
        intracellular_stim.amp = fiber.search(Config.SIM, 'intracellular_stim', 'amp')
        self.istim = intracellular_stim
        return self

    def initialize_extracellular(self, fiber):
        pass

    def update_extracellular(self, fiber, e_stims):
        pass