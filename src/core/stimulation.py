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
        if fiber.fiber_type == 2:
            for sec in fiber.node:
                sec(0.5).e_extracellular = 0
            for sec in fiber.MYSA:
                sec(0.5).e_extracellular = 0
            for sec in fiber.FLUT:
                sec(0.5).e_extracellular = 0
            for sec in fiber.STIN:
                sec(0.5).e_extracellular = 0
        else:
            for sec in fiber.sec:
                sec(0.5).e_extracellular = 0
        return 

    def update_extracellular(self, fiber, e_stims):
        if fiber.fiber_type == 2:
            node_stim, FLUT_stim, MYSA_stim, STIN_stim = [], [], [], []
            for ind in range(1, len(e_stims) + 1):
                if ind % 11 == 1:
                    node_stim.append(e_stims[ind - 1])
                elif ind % 11 == 2 or ind % 11 == 0:
                    MYSA_stim.append(e_stims[ind - 1])
                elif ind % 11 == 3 or ind % 11 == 10:
                    FLUT_stim.append(e_stims[ind - 1])
                else:
                    STIN_stim.append(e_stims[ind - 1])

            for x, sec in enumerate(fiber.node):
                sec(0.5).e_extracellular = node_stim[x]
            for x, sec in enumerate(fiber.MYSA):
                sec(0.5).e_extracellular = MYSA_stim[x]
            for x, sec in enumerate(fiber.FLUT):
                sec(0.5).e_extracellular = FLUT_stim[x]
            for x, sec in enumerate(fiber.STIN):
                sec(0.5).e_extracellular = STIN_stim[x]
        else:
            for x, sec in enumerate(fiber.sec):
                sec(0.5).e_extracellular = e_stims[x]
