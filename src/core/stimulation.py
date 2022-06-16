"""
Extracellular Stimulation class

MISSING: error checks from original .hoc files
"""
from neuron import h
import pickle
h.load_file('stdrun.hoc')

class ExtracellularStimulation():
    def __init__(self):
        self.VeSpace_data = []
        self.VeTime_data =[]

    def load_space(self, space_file):
        """
        :param fname: file name containing Extracellular Stim Space data
        :return: Ve(x) -- vector of potentials from FEM
        """
        spacefile = open(space_file, 'r')
        self.axontotal = int(spacefile.readline())
        file_lines = spacefile.read().splitlines()
        VeSpace_data = [float(i)*1000 for i in file_lines]  # Need to convert to V -> mV
        spacefile.close()

        if len(VeSpace_data) != self.axontotal:
            raise Exception("Need axontotal from VeSpace file to match axontotal used in Python")
        self.VeSpace_data = VeSpace_data
        return self


    def load_time(self, time_file):
        """
        :param fname: file name containing Extracellular Stim Time data
        :return: I(t) -- vector of amplitudes at each time step of the FEM
        """
        timefile = open(time_file, 'r')
        self.dt = float(timefile.readline().strip())           # time step
        self.tstop = int(timefile.readline().strip())          # stop time
        file_lines = timefile.read().splitlines()
        VeTime_data = [float(i) for i in file_lines]
        timefile.close()

        self.VeTime_data = VeTime_data
        return self
