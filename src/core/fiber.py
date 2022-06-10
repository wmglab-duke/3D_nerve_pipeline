"""
Fiber class: fiber object with relevant attributes
"""
import json

from neuron import h

from src.utils import (Config, Configurable, Exceptionable, MyelinationMode, SetupMode, WriteMode)
from src.core.extracellularStim import ExtracellularStimulation

import math
import numpy as np
import os
import bokeh.plotting as plt
from bokeh.layouts import row
import pickle

Section = h.Section
SectionList = h.SectionList
ion_style = h.ion_style

h.load_file('stdrun.hoc')


""" For an unmyelinated fiber, return a list of variable density, 
such that the center of the fiber has greater density than the edges"""
def dz_calculator() -> list:
    return []


class Fiber(Exceptionable, Configurable):
    def __init__(self):
        exceptions_file = os.path.join('config', 'system', 'exceptions.json')

        with open(exceptions_file, "r") as handle:
            exceptions_config: dict = json.load(handle)

        Configurable.__init__(self)
        Exceptionable.__init__(self, SetupMode.OLD, exceptions_config)

        self.xyz = tuple()
        self.index = None
        self.diameter = None
        self.min = None
        self.max = None
        self.offset = None
        self.seed = None
        self.myelination = None
        self.delta_z = None
        self.paranodal_length_2 = None
        self.inter_length = None
        self.passive_end_nodes = None
        self.node = []
        self.MYSA = []
        self.FLUT = []
        self.STIN = []
        self.last_run = bool
        self.v_init = None
        self.v_init_c_fiber = None
        return

    def inherit(self, xyz: tuple, index: int):
        """
        Inherit known properties of the fiber based on sim config
        """
        self.index = index
        self.xyz = xyz
        self.diameter = self.search(Config.SIM, 'fibers', 'z_parameters', 'diameter')
        self.min = self.search(Config.SIM, 'fibers', 'z_parameters', 'min')
        self.max = self.search(Config.SIM, 'fibers', 'z_parameters', 'max')
        self.offset = self.search(Config.SIM, 'fibers', 'z_parameters', 'offset')
        self.seed = self.search(Config.SIM, 'fibers', 'z_parameters', 'seed')
        self.fiber_mode = self.search(Config.SIM, 'fibers', 'mode')
        self.myelination = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'myelinated')
        self.temperature = self.search(Config.MODEL, 'temperature')
        return self

    def generate(self, fiber_path):
        """
        Build fiber sections based on fiber type
        Reads in geometric properties from JSON files
        Makes calls to construct fiber sections
        """

        if self.fiber_mode != 'MRG_DISCRETE' and self.fiber_mode != 'MRG_INTERPOLATION':
            self.node_channels = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'node_channels')
            self.delta_z = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'delta_zs')
            self.fiber_type = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'fiber_type')
            self.passive_end_nodes = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'passive_end_nodes')
            self.channels_type = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'channels_type')
            self.neuron_flag = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'neuron_flag')


        elif self.fiber_mode == 'MRG_DISCRETE':
            diameters, my_delta_zs, paranodal_length_2s = (
                self.search(Config.FIBER_Z, MyelinationMode.parameters.value, self.fiber_mode, key)
                for key in ('diameters', 'delta_zs', 'paranodal_length_2s')
            )
            diameter_index = diameters.index(self.diameter)
            self.node_channels = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'node_channels')
            delta_z = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'delta_zs')[diameter_index]
            node_length = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'node_length')
            paranodal_length_1 = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                             'paranodal_length_1')
            paranodal_length_2 = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                             'paranodal_length_2s')[diameter_index]
            self.neuron_flag = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'neuron_flag')
            self.fiber_type = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'fiber_type')
            self.inter_length = eval(self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'inter_length'))
            self.fiber_type = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'fiber_type')
            self.passive_end_nodes = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'passive_end_nodes')
            self.delta_z, self.node_length, self.paranodal_length_1, self.paranodal_length_2 = delta_z, node_length, \
                                                                                               paranodal_length_1, paranodal_length_2

            if self.diameter == 1:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = None, 0.8, 0.7, 0.7, 0.8, 15
            elif self.diameter == 2:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = None, 1.6, 1.4, 1.4, 1.6, 30
            elif self.diameter == 5.7:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = 0.605, 3.4, 1.9, 1.9, 3.4, 80
            elif self.diameter == 7.3:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = 0.630, 4.6, 2.4, 2.4, 4.6, 100
            elif self.diameter == 8.7:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = 0.661, 5.8, 2.8, 2.8, 5.8, 110
            elif self.diameter == 10:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = 0.690, 6.9, 3.3, 3.3, 6.9, 120
            elif self.diameter == 11.5:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = 0.700, 8.1, 3.7, 3.7, 8.1, 130
            elif self.diameter == 12.8:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = 0.719, 9.2, 4.2, 4.2, 9.2, 135
            elif self.diameter == 14:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = 0.739, 10.4, 4.7, 4.7, 10.4, 140
            elif self.diameter == 15:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = 0.767, 11.5, 5.0, 5.0, 11.5, 145
            elif self.diameter == 16:
                self.g, self.axonD, self.nodeD, self.paraD1, self.paraD2, self.nl = 0.791, 12.7, 5.5, 5.5, 12.7, 150

        elif self.fiber_mode == 'MRG_INTERPOLATION':
            self.node_length = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'node_length')
            self.paranodal_length_1 = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                             'paranodal_length_1')
            self.paranodal_length_2 = eval(self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                             'paranodal_length_2'))

            self.inter_length = eval(self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                             'inter_length'))
            self.fiber_type = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'fiber_type')
            self.passive_end_nodes = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                            'passive_end_nodes')
            if self.diameter >= 5.643:
                self.delta_z = eval(self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                           'delta_z', 'diameter_greater_or_equal_5.643um'))
            else:
                self.delta_z = eval(self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                           'delta_z', 'diameter_less_5.643um'))

            self.nl = -0.4749 * self.diameter ** 2 + 16.85 * self.diameter - 0.7648
            self.nodeD = 0.01093 * self.diameter ** 2 + 0.1008 * self.diameter + 1.099
            self.paraD1 = self.nodeD
            self.paraD2 = 0.02361 * self.diameter ** 2 + 0.3673 * self.diameter + 0.7122
            self.axonD = self.paraD2


        fiber_ve = np.loadtxt(fiber_path)
        n_fiber_coords = int(fiber_ve[0])

        if self.neuron_flag == 2:
            self.axonnodes = int(1 + (n_fiber_coords - 1) / 11)
        elif self.neuron_flag == 3:
            self.axonnodes = int(n_fiber_coords)
            self.len = self.delta_z*self.axonnodes

        if self.myelination:
            self.createMyelinatedFiber(self.node_channels, self.axonnodes, self.diameter, self.temperature, self.axonD,
                                       self.nodeD, self.paraD1, self.paraD2, self.delta_z, self.paranodal_length_2,
                                       self.nl, self.passive_end_nodes)
        elif not self.myelination:
            self.createUnmyelinatedFiber(self.diameter, self.len, c_fiber_model_type=self.channels_type, celsius=self.temperature,
                                         delta_z=self.delta_z, passive_end_nodes=self.passive_end_nodes)

        return self

    def createMyelinatedFiber(self, node_channels, axonnodes, fiberD, celsius, axonD, nodeD, paraD1, paraD2, deltaz,
                              paralength2, nl, passive_end_nodes):
        # Electrical parameters
        rhoa = 0.7e6        # [ohm-um]
        mycm = 0.1          # lamella membrane; [uF/cm2]
        mygm = 0.001        # lamella membrane; [S/cm2]
        rhoe = 1000         # resistivity of extracellular medium; [ohm-cm]

        if node_channels == 0:
            e_pas_Vrest = -80
        elif node_channels == 1:
            e_pas_Vrest = -57

        # Geometrical parameters [um]
        paranodes1 = 2*(axonnodes-1)    # Number of MYSA paranodes
        paranodes2 = 2*(axonnodes-1)    # Number of FLUT paranodes
        axoninter = 6*(axonnodes-1)     # Number of STIN internodes

        nodelength = 1.0                # Length of nodes of Ranvier [um]
        paralength1 = 3                 # Length of MYSA [um]
        space_p1 = 0.002                # Thickness of periaxonal space in MYSA sections [um]
        space_p2 = 0.004                # Thickness of periaxonal space in FLUT sections [um]
        space_i = 0.004                 # Thickness of periaxonal space in STIN sections [um]

        Rpn0 = (rhoa * .01) / (math.pi * ((((nodeD / 2) + space_p1)**2) - ((nodeD / 2)**2)))
        Rpn1 = (rhoa * .01) / (math.pi * ((((paraD1 / 2) + space_p1)**2) - ((paraD1 / 2)**2)))
        Rpn2 = (rhoa * .01) / (math.pi * ((((paraD2 / 2) + space_p2)**2) - ((paraD2 / 2)**2)))
        Rpx = (rhoa * .01) / (math.pi * ((((axonD / 2) + space_i)**2) - ((axonD / 2)**2)))
        interlength = (deltaz - nodelength - (2 * paralength1) - (2 * paralength2)) / 6

        # Create the axon sections
        for i in range(0, axonnodes):
            new_node = node_creator(i, nodeD, nodelength, rhoa, mycm, mygm, rhoe, passive_end_nodes, axonnodes,
                                    node_channels, nl, Rpn0, celsius)
            self.node.append(new_node.node)
        for i in range(0, paranodes1):
            new_MYSA = MYSA_creator(i, fiberD, paralength1, rhoa, paraD1, e_pas_Vrest, Rpn1, mycm, mygm, nl)
            self.MYSA.append(new_MYSA.obj)
        for i in range(0, paranodes2):
            new_FLUT = FLUT_creator(i, fiberD, paralength2, rhoa, paraD2, e_pas_Vrest, Rpn2, mycm, mygm, nl)
            self.FLUT.append(new_FLUT.obj)
        for i in range(0, axoninter):
            new_STIN = STIN_creator(i, fiberD, interlength, rhoa, axonD, e_pas_Vrest, Rpx, mycm, mygm, nl)
            self.STIN.append(new_STIN.obj)

        # Connect the axon sections
        for i in range(0, axonnodes-1):
            self.MYSA[2*i].connect(self.node[i])
            self.FLUT[2*i].connect(self.MYSA[2*i])
            self.STIN[6*i].connect(self.FLUT[2*i])
            self.STIN[6*i+1].connect(self.STIN[6*i])
            self.STIN[6*i+2].connect(self.STIN[6*i+1])
            self.STIN[6*i+3].connect(self.STIN[6*i+2])
            self.STIN[6*i+4].connect(self.STIN[6*i+3])
            self.STIN[6*i+5].connect(self.STIN[6*i+4])
            self.FLUT[2*i+1].connect(self.STIN[6*i+5])
            self.MYSA[2*i+1].connect(self.FLUT[2*i+1])
            self.node[i+1].connect(self.MYSA[2*i+1])
        return self

    def createUnmyelinatedFiber(self, fiberD=6, len=21, c_fiber_model_type=1, celsius=37, delta_z=50/6, insert97na=0,
                     conductances97=0, passive_end_nodes=0):
        """
        Create a list of Neuron Sections for an unmyelinated fiber
        """
        nsegments = int(len/delta_z)

        self.sec = []
        for i in range(0, nsegments):
            node = Section(name='node ' + str(i))
            self.sec.append(node)

            node.diam = fiberD
            node.nseg = 1
            node.L = delta_z
            if passive_end_nodes:
                if i == 0 or i == nsegments-1:
                    node.insert('pas')
                    node.g_pas = 0.0001
                    if c_fiber_model_type == 1:
                        node.e_pas = -60       # Sundt model equilibrium potential
                    elif c_fiber_model_type == 2:
                        node.e_pas = -55       # Tigerholm model equilibrium potential
                    elif c_fiber_model_type == 3:
                        node.e_pas = -70       # Rattay model equilibrium potential
                    elif c_fiber_model_type == 4:
                        node.e_pas = -48       # Schild model equilibrium potential
                    else:
                        node.e_pas = -70

                node.insert('extracellular')
                node.xc[0] = 0  # short circuit
                node.xg[0] = 1e10  # short circuit

                node.Ra = 1e10

            else:
                if c_fiber_model_type == 1:     # Sundt model
                    node.insert('nahh')
                    node.gnabar_nahh = 0.04
                    node.mshift_nahh = -6          # NaV1.7/1.8 channelshift
                    node.hshift_nahh = 6           # NaV1.7/1.8 channelshift

                    node.insert('borgkdr')         # insert delayed rectified K channels
                    node.gkdrbar_borgkdr = 0.04    # density of K channels
                    node.ek = -90                  # K equilibrium potential

                    node.insert('pas')             # insert leak channels
                    node.g_pas = 1/10000           # set Rm = 10000 ohms-cm2
                    node.Ra = 100                  # intracellular resistance
                    node.v = -60
                    node.e_pas = node.v + (node.ina + node.ik)/node.g_pas      # calculate leak equilibrium potential
                elif c_fiber_model_type == 2:   # Tigerholm model
                    node.insert('ks')
                    node.gbar_ks = 0.0069733
                    node.insert('kf')
                    node.gbar_kf = 0.012756
                    node.insert('h')
                    node.gbar_h = 0.0025377
                    node.insert('nattxs')
                    node.gbar_nattxs = 0.10664
                    node.insert('nav1p8')
                    node.gbar_nav1p8 = 0.24271
                    node.insert('nav1p9')
                    node.gbar_nav1p9 = 9.4779e-05
                    node.insert('nakpump')
                    node.smalla_nakpump = -0.0047891
                    node.insert('kdrTiger')
                    node.gbar_kdrTiger = 0.018002
                    node.insert('kna')
                    node.gbar_kna = 0.00042
                    node.insert('naoiTiger')
                    # node.theta_naoi	= 0.029
                    node.insert('koiTiger')
                    # node.theta_koi	= 0.029

                    node.insert('leak')
                    node.insert('extrapump')

                    node.Ra = 35.5
                    node.cm = 1

                    node.celsiusT_ks = celsius
                    node.celsiusT_kf = celsius
                    node.celsiusT_h = celsius
                    node.celsiusT_nattxs = celsius
                    node.celsiusT_nav1p8 = celsius
                    node.celsiusT_nav1p9 = celsius
                    node.celsiusT_nakpump = celsius
                    node.celsiusT_kdrTiger = celsius
                    node.v = -55
                elif c_fiber_model_type == 3:  # Rattay and Aberham model -- adjusted for a resting potential of -70mV
                    node.insert('RattayAberham')
                    node.Ra = 100      # required for propagation; less than 100 does not propagate
                    node.cm = 1
                    node.v = -70
                    node.ena = 45
                    node.ek = -82
                elif c_fiber_model_type == 4:   # Schild model
                    node.R = 8314
                    node.F = 96500
                    node.insert('leakSchild')						# All mechanisms from Schild 1994 inserted into model
                    node.insert('kd')
                    node.insert('ka')
                    node.insert('can')
                    node.insert('cat')
                    node.insert('kds')
                    node.insert('kca')
                    node.insert('caextscale')
                    node.insert('caintscale')
                    node.insert('CaPump')
                    node.insert('NaCaPump')
                    node.insert('NaKpumpSchild')
                    if insert97na:
                        node.insert('naf97mean')
                        node.insert('nas97mean')
                    else:
                        node.insert('naf')
                        node.insert('nas')

                    # Ionic concentrations
                    node.cao0_ca_ion = 2.0                                               # [mM] Initial Cao Concentration
                    node.cai0_ca_ion = 0.000117                                          # [mM] Initial Cai Concentrations
                    node.ko = 5.4                                                        # [mM] External K Concentration
                    node.ki = 145.0                                                      # [mM] Internal K Concentration
                    node.kstyle = ion_style("k_ion", 1, 2, 0, 0, 0)                      # Allows ek to be calculated manually
                    node.ek = ((node.R*(celsius+273.15))/node.F)*np.log10(node.ko/node.ki)           # Manual Calculation of ek in order to use Schild F and R values

                    node.nao = 154                                                       # [mM] External Na Concentration
                    node.nai = 8.9                                                       # [mM] Internal Na Concentration
                    node.nastyle = ion_style("na_ion", 1, 2, 0, 0, 0)                    # Allows ena to be calculated manually
                    node.ena = ((node.R*(celsius+273.15))/node.F)*np.log10(node.nao/node.nai)        # Manual Calculation of ena in order to use Schild F and R values
                    if conductances97:
                        node.gbar_naf97mean = 0.022434928                                # [S/cm^2] This block sets the conductance to the conductances in Schild 1997
                        node.gbar_nas97mean = 0.022434928
                        node.gbar_kd = 0.001956534
                        node.gbar_ka = 0.001304356
                        node.gbar_kds = 0.000782614
                        node.gbar_kca = 0.000913049
                        node.gbar_can = 0.000521743
                        node.gbar_cat = 0.00018261
                        node.gbna_leak = 1.8261E-05
                    node.Ra = 100
                    node.cm = 1.326291192
                    node.v = -48

                self.v_init_c_fiber = node.v

                node.insert('extracellular')
                node.xc[0] = 0  # short circuit
                node.xg[0] = 1e10  # short circuit

        for i in range(0, nsegments-1):
            self.sec[i+1].connect(self.sec[i])

        return self

    def balance(self):
        Vrest = -55
        for s in self.sec:
            if (-(s.ina_nattxs + s.ina_nav1p9 + s.ina_nav1p8 + s.ina_h + s.ina_nakpump)/(Vrest - s.ena)) < 0:
                s.pumpina_extrapump = -(s.ina_nattxs + s.ina_nav1p9 + s.ina_nav1p8 + s.ina_h + s.ina_nakpump)
            else:
                s.gnaleak_leak = -(s.ina_nattxs + s.ina_nav1p9 + s.ina_nav1p8 + s.ina_h + s.ina_nakpump) / (Vrest - s.ena)

            if (-(s.ik_ks + s.ik_kf + s.ik_h + s.ik_kdrTiger + s.ik_nakpump + s.ik_kna) / (Vrest - s.ek)) < 0:
                s.pumpik_extrapump = -(s.ik_ks + s.ik_kf + s.ik_h + s.ik_kdrTiger + s.ik_nakpump + s.ik_kna)
            else:
                s.gkleak_leak = -(s.ik_ks + s.ik_kf + s.ik_h + s.ik_kdrTiger + s.ik_nakpump + s.ik_kna) / (Vrest - s.ek)

    def write(self, mode: WriteMode, path: str):
        """
        :param mode:
        :param path:
        :return:
        """

        diams = []
        with open(os.path.join(path, WriteMode.file_endings.value[mode.value]), 'w') as f:
            for row in [len(self.z)] + list(self.z):
                if not isinstance(row, int):
                    for el in row:
                        f.write(str(el) + ' ')
                else:
                    f.write(str(row) + ' ')
                f.write("\n")
        return self

    def run(self, stimamp, fiber_path, waveform_path, n_tsteps, plot=False):
        if plot:
            ap_detect_location = self.search(Config.SIM, 'protocol', 'threshold', 'ap_detect_location')
            node_index = int((self.axonnodes - 1) * self.delta_z * ap_detect_location / self.delta_z)
            plot_node = plt.figure(title='node[' + str(node_index) + '] at stimamp = ' + str(stimamp),
                                   x_axis_label='t (ms)', y_axis_label='v (mV)')
            v_node = h.Vector().record(self.sec[node_index](0.5)._ref_v)
            t = h.Vector().record(h._ref_t)

        for sec in self.sec:
            if sec.v != -55:
                print('sec.v wrong prior to initialize')
                break

        h.finitialize(self.v_init)

        for sec in self.sec:
            if sec.v != -55:
                print('sec.v wrong after initialize')
                break

        if self.fiber_type == 3:
            if self.channels_type == 2 and self.passive_end_nodes == 1:
                raise Exception("Program cannot balance Tigerholm for passive_end_nodes=1, must be 0.")
            elif self.channels_type == 2 and self.passive_end_nodes == 0:
                self.balance()

        for sec in self.sec:
            if sec.v != -55:
                print('sec.v wrong after balance')
                break

        stimulation_vectors = ExtracellularStimulation()
        stimulation_vectors \
            .load_space(fiber_path) \
            .load_time(waveform_path)
        t_step = stimulation_vectors.dt
        tstop = stimulation_vectors.tstop

        # Determine which gating parameters to save
        space_gating: bool = self.search(Config.SIM, 'saving', 'space', 'gating')
        time_gating: bool = self.search(Config.SIM, 'saving', 'time', 'gating')
        runtimes: bool = self.search(Config.SIM, 'saving', 'runtimes')

        # Determine protocol
        protocol_mode = self.search(Config.SIM, 'protocol', 'mode')
        t_initSS = self.search(Config.SIM, 'protocol', 'initSS')
        dt_initSS = self.search(Config.SIM, 'protocol', 'dt_initSS')

        if self.fiber_type == 2:
            for sec in self.node:
                sec(0.5).e_extracellular = 0
            for sec in self.MYSA:
                sec(0.5).e_extracellular = 0
            for sec in self.FLUT:
                sec(0.5).e_extracellular = 0
            for sec in self.STIN:
                sec(0.5).e_extracellular = 0
        else:
            for sec in self.sec:
                sec(0.5).e_extracellular = 0

        h.t = t_initSS      # Start before t=0
        h.dt = dt_initSS    # Large dt
        while (h.t <= -dt_initSS):
            h.fadvance()

        h.dt = t_step
        h.t = 0
        h.fcurrent()
        h.frecord_init()
        h.celsius = self.temperature

        # Set up APcount
        apc = []
        if self.fiber_type == 2:
            for i, node in enumerate(self.node):
                apc.append(h.APCount(node(0.5)))
                apc[i].thresh = self.search(Config.SIM, "protocol", "threshold", "value")
        else:
            for i, node in enumerate(self.sec):
                apc.append(h.APCount(node(0.5)))
                apc[i].thresh = self.search(Config.SIM, "protocol", "threshold", "value")

        # Begin time loop
        for i in range(0, n_tsteps):
            if i%(int(n_tsteps/5)) == 0:
               print(str(i) + '/' + str(n_tsteps) + " time steps")
            if i*t_step > tstop:
                break
            amp = stimulation_vectors.VeTime_data[i]
            scaled_stim = [stimamp*amp*j for j in stimulation_vectors.VeSpace_data]
            if self.fiber_type == 2:
                node_stim, FLUT_stim, MYSA_stim, STIN_stim = [], [], [], []
                for ind in range(1, len(scaled_stim)+1):
                    if ind%11 == 1:
                        node_stim.append(scaled_stim[ind-1])
                    elif ind%11 == 2 or ind%11 == 0:
                        MYSA_stim.append(scaled_stim[ind-1])
                    elif ind%11 == 3 or ind%11 == 10:
                        FLUT_stim.append(scaled_stim[ind-1])
                    else:
                        STIN_stim.append(scaled_stim[ind-1])

                for x, sec in enumerate(self.node):
                    sec(0.5).e_extracellular = node_stim[x]
                for x, sec in enumerate(self.MYSA):
                    sec(0.5).e_extracellular = MYSA_stim[x]
                for x, sec in enumerate(self.FLUT):
                    sec(0.5).e_extracellular = FLUT_stim[x]
                for x, sec in enumerate(self.STIN):
                    sec(0.5).e_extracellular = STIN_stim[x]
            else:
                for x, sec in enumerate(self.sec):
                    sec(0.5).e_extracellular = scaled_stim[x]
            h.fadvance()

        if plot:
            plot_node.line(t, list(v_node), line_width=2)
            plt.show(plot_node)

            outfile = open('model_tiger_data', 'wb')
            pickle.dump([t, list(v_node)], outfile)
            outfile.close()


        ap_detect_location = self.search(Config.SIM, 'protocol', 'threshold', 'ap_detect_location')
        node_index = int((self.axonnodes-1)*self.delta_z*ap_detect_location/self.delta_z)
        if apc[node_index].n >= 1:
            self.last_run = True
        else:
            self.last_run = False
        return self


class GeometryObject():
    def __init__(self, fiberD, fiberDtoAxonD=0, axonDtoNL=0, nodelength=1, MYSAlength=3):
        self.fiberDtoAxonD, self.axonDtoNL, self.nodelength, self.MYSAlength = fiberDtoAxonD, axonDtoNL, nodelength, MYSAlength
        self.FLUTlength = -0.171 * (fiberD**2) + 6.48 * fiberD - 0.935
        if fiberDtoAxonD == 0:
            self.axonD = 0.553 * fiberD - 0.024
        elif fiberDtoAxonD == 1:
            self.axonD = 0.688 * fiberD - 0.337
        elif fiberDtoAxonD == 2:
            self.axonD = 0.0156 * (fiberD**2) + 0.392 * fiberD + 0.188
        elif fiberDtoAxonD == 3:
            self.axonD = 0.684 * fiberD + 0.0821
        elif fiberDtoAxonD == 4:
            self.axonD = 0.621 * fiberD - 0.121

        self.nodeD = 0.321 * self.axonD + 0.37

        if axonDtoNL == 0:
            self.nl = int(17.4 * self.axonD - 1.74)
        elif axonDtoNL == 1:
            self.nl = int(-1.17 * (self.axonD**2) + 24.9 * self.axonD + 17.7)

        self.MYSAD = self.nodeD
        self.FLUTD = self.axonD

        self.interlength = ((-3.22 * fiberD**2 + 148 * fiberD + -128) - self.nl - 2*self.MYSAlength - 2*self.FLUTlength)/6

class MYSA_creator():
    def __init__(self, i, fiberD, paralength1, rhoa, paraD1, e_pas_Vrest, Rpn1, mycm, mygm, nl):
        self.obj = Section(name='MYSA ' + str(i))
        self.obj.nseg = 1
        self.obj.diam = fiberD
        self.obj.L = paralength1
        self.obj.Ra = rhoa * (1 / (paraD1 / fiberD) ** 2) / 10000
        self.obj.cm = 2 * paraD1 / fiberD
        self.obj.insert('pas')
        self.obj.g_pas = 0.001 * paraD1 / fiberD
        self.obj.e_pas = e_pas_Vrest

        self.obj.insert('extracellular')
        self.obj.xraxial[0] = Rpn1
        self.obj.xc[0] = mycm / (nl * 2)  # short circuit
        self.obj.xg[0] = mygm / (nl * 2)  # short circuit
        # for ind in range(0, 2):
        #     self.obj.xraxial[ind] = Rpn1
        #     self.obj.xc[ind] = mycm / (nl * 2)  # short circuit
        #     self.obj.xg[ind] = mygm / (nl * 2)  # short circuit
        return

class FLUT_creator():
    def __init__(self, i, fiberD, paralength2, rhoa, paraD2, e_pas_Vrest, Rpn2, mycm, mygm, nl):
        self.obj = Section(name='FLUT ' + str(i))
        self.obj.nseg = 1
        self.obj.diam = fiberD
        self.obj.L = paralength2
        self.obj.Ra = rhoa * (1 / (paraD2 / fiberD) ** 2) / 10000
        self.obj.cm = 2 * paraD2 / fiberD
        self.obj.insert('pas')
        self.obj.g_pas = 0.0001 * paraD2 / fiberD
        self.obj.e_pas = e_pas_Vrest

        self.obj.insert('extracellular')
        self.obj.xraxial[0] = Rpn2
        self.obj.xc[0] = mycm / (nl * 2)  # short circuit
        self.obj.xg[0] = mygm / (nl * 2)  # short circuit

        # for ind in range(0, 2):
        #     self.obj.xraxial[ind] = Rpn2
        #     self.obj.xc[ind] = mycm / (nl * 2)  # short circuit
        #     self.obj.xg[ind] = mygm / (nl * 2)  # short circuit

        return

class STIN_creator():
    def __init__(self, i, fiberD, interlength, rhoa, axonD, e_pas_Vrest, Rpx, mycm, mygm, nl):
        self.obj = Section(name='STIN ' + str(i))
        self.obj.nseg = 1
        self.obj.diam = fiberD
        self.obj.L = interlength
        self.obj.Ra = rhoa * (1 / (axonD / fiberD) ** 2) / 10000
        self.obj.cm = 2 * axonD / fiberD
        self.obj.insert('pas')
        self.obj.g_pas = 0.0001 * axonD / fiberD
        self.obj.e_pas = e_pas_Vrest

        self.obj.insert('extracellular')
        self.obj.xraxial[0] = Rpx
        self.obj.xc[0] = mycm / (nl * 2)  # short circuit
        self.obj.xg[0] = mygm / (nl * 2)  # short circuit

        # for ind in range(0, 2):
        #     self.obj.xraxial[ind] = Rpx
        #     self.obj.xc[ind] = mycm / (nl * 2)  # short circuit
        #     self.obj.xg[ind] = mygm / (nl * 2)  # short circuit

        return

class node_creator():
    def __init__(self, index, nodeD, nodelength, rhoa, mycm, mygm, rhoe, passive, axonnodes, node_channels, nl, Rpn0, celsius):
        self.node = Section(name='node ' + str(index))
        self.node.nseg = 1
        self.node.diam = nodeD
        self.node.L = nodelength
        self.node.Ra = rhoa/10000

        if passive and (index == 0 or index == axonnodes - 1):
            self.node.cm = 2
            self.node.insert('pas')
            self.node.g_pas = 0.0001
            self.node.e_pas = -70
            self.node.insert('extracellular')
            self.node.xc[0] = mycm / (nl * 2)  # short circuit
            self.node.xg[0] = mygm / (nl * 2)  # short circuit
            # for ind in range(0, 2):
            #     self.node.xc[ind] = mycm / (nl * 2)  # short circuit
            #     self.node.xg[ind] = mygm / (nl * 2)  # short circuit

        else:
            if node_channels == 0:
                self.node.cm = 2
                self.node.insert('axnode_myel')

            elif node_channels == 1:
                self.node.cm = 1.149452367   # [uF/cm^2] specific membrane capacitance (Schild 1994, A-type)
                F = 96500                    # [C/mole] Faraday'node Constant from Schild 1994
                R = 8314                     # [J/(kg*mole*K)] Gas Constant from Schild 1994

                # Based on Schild 1994 ion channels
                self.node.insert('leakSchild')
                self.node.insert('naf')
                self.node.insert('nas')
                self.node.insert('kd')
                self.node.insert('ka')
                self.node.insert('can')
                self.node.insert('cat')
                self.node.insert('kds')
                self.node.insert('kca')

                self.node.insert('caextscale')
                self.node.insert('caintscale')
                self.node.insert('CaPump')
                self.node.insert('NaCaPump')
                self.node.insert('NakpumpSchild')

                self.node.L_caintscale = self.node.L
                self.node.nseg_caintscale = self.node.nseg
                self.node.L_caextscale = self.node.L
                self.node.nseg_caextscale = self.node.nseg

                # Ionic concentrations
                self.node.cao0_ca_ion = 2.0  # [mM] Initial Cao Concentration
                self.node.cai0_ca_ion = 0.000117  # [mM] Initial Cai Concentrations
                self.node.ko = 5.4  # [mM] External K Concentration
                self.node.ki = 145.0  # [mM] Internal K Concentration
                kstyle = ion_style("k_ion", 1, 2, 0, 0, 0)  # Allows ek to be calculated manually
                self.node.ek = ((R * (celsius + 273.15)) / F) * np.log10(
                    self.node.ko / self.node.ki)  # Manual Calculation of ek in order to use Schild F and R values

                self.node.nao = 154  # [mM] External Na Concentration
                self.node.nai = 8.9  # [mM] Internal Na Concentration
                self.nastyle = ion_style("na_ion", 1, 2, 0, 0, 0)  # Allows ena to be calculated manually
                self.node.ena = ((R * (celsius + 273.15)) / F) * np.log10(
                    self.node.nao / self.node.nai)  # Manual Calculation of ena in order to use Schild F and R values

                self.node.gbar_naf = 3 // 0.072503919  # NOTE: Does not conduct with original Schild 1994 value of 0.072503919; does not conduct at 0.5, but does at 1; increased to MRG value of 3
                self.node.shiftnaf_naf = 0  # [mV]
                self.node.gbar_nas = 3.53678E-07
                self.node.shiftnas_nas = 0
                self.node.gbar_kd = 0.000194523
                self.node.shiftkd_kd = 0
                self.node.gbar_ka = 0.001237872
                self.node.shiftka_ka = 0
                self.node.gbar_kds = 0.000353678
                self.node.shiftkds_kds = 0
                self.node.gbar_kca = 0.00022989
                self.node.gbar_can = 3.53678E-05
                self.node.shiftcan_can = 0
                self.node.gbar_cat = 1.23787E-05
                self.node.shiftcan_cat = 0
                self.node.gbna_leak = 1.14945E-05

            self.node.insert('extracellular')
            self.node.xraxial[0] = Rpn0
            self.node.xc[0] = 0  # short circuit
            self.node.xg[0] = 1e10  # short circuit
            # for ind in range(0, 2):
            #     self.node.xraxial[ind] = Rpn0
            #     self.node.xc[ind] = 0         # short circuit
            #     self.node.xg[ind] = 1e10      # short circuit
        return

