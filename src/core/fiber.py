from neuron import h
from src.utils.enums import (Config, MyelinationMode, NeuronRunMode, TerminationCriteriaMode, SearchAmplitudeIncrementMode)
from src.utils.configurable import Configurable
from src.utils.saveable import Saveable

import math
import numpy as np
import time

Section = h.Section
SectionList = h.SectionList
ion_style = h.ion_style

h.load_file('stdrun.hoc')

""" For an unmyelinated fiber, return a list of variable density, 
such that the center of the fiber has greater density than the edges"""
def dz_calculator() -> list:
    return []


class Fiber(Configurable, Saveable):
    def __init__(self):
        Configurable.__init__(self)

        self.diameter = None
        self.min = None
        self.max = None
        self.offset = None
        self.seed = None
        self.fiber_mode = None
        self.myelination = None
        self.temperature = None
        self.axonnodes = None
        self.delta_z = None
        self.passive_end_nodes = None
        self.node = []
        self.MYSA = []
        self.FLUT = []
        self.STIN = []
        self.sec = []
        self.n_aps = None
        self.inner_ind = None
        self.fiber_ind = None
        self.v_init = None
        return

    def inherit(self):
        """
        Inherit known properties of the fiber based on sim config
        """
        self.diameter = self.search(Config.SIM, 'fibers', 'z_parameters', 'diameter')
        self.min = self.search(Config.SIM, 'fibers', 'z_parameters', 'min')
        self.max = self.search(Config.SIM, 'fibers', 'z_parameters', 'max')
        self.offset = self.search(Config.SIM, 'fibers', 'z_parameters', 'offset')
        self.seed = self.search(Config.SIM, 'fibers', 'z_parameters', 'seed')
        self.fiber_mode = self.search(Config.SIM, 'fibers', 'mode')
        self.myelination = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'myelinated')
        self.temperature = self.search(Config.MODEL, 'temperature')
        return self

    def generate(self, n_fiber_coords):
        """
        Build fiber sections based on fiber type
        Reads in geometric properties from JSON files
        Makes calls to construct fiber sections
        """

        if self.fiber_mode != 'MRG_DISCRETE' and self.fiber_mode != 'MRG_INTERPOLATION':
            fiber_type = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'fiber_type')
            neuron_flag = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'neuron_flag')
            node_channels = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'node_channels')
            self.delta_z = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'delta_zs')
            self.passive_end_nodes = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'passive_end_nodes')
            channels_type = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'channels_type')

        elif self.fiber_mode == 'MRG_DISCRETE':
            diameters, my_delta_zs, paranodal_length_2s, gs, axonDs, nodeDs, paraD1s, paraD2s, nls = (
                self.search(Config.FIBER_Z, MyelinationMode.parameters.value, self.fiber_mode, key)
                for key in ('diameters', 'delta_zs', 'paranodal_length_2s', "gs", "axonDs",
                            "nodeDs", "paraD1s", "paraD2s", "nls")
            )
            diameter_index = diameters.index(self.diameter)
            neuron_flag = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'neuron_flag')
            node_channels = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'node_channels')
            self.delta_z = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'delta_zs')[diameter_index]
            paranodal_length_2 = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                             'paranodal_length_2s')[diameter_index]
            g = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'gs')[diameter_index]
            axonD = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'axonDs')[diameter_index]
            nodeD = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'axonDs')[diameter_index]
            paraD1 = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'paraD1s')[diameter_index]
            paraD2 = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'paraD2s')[diameter_index]
            nl = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'nls')[diameter_index]
            self.passive_end_nodes = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'passive_end_nodes')
            fiber_type = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'fiber_type')

        elif self.fiber_mode == 'MRG_INTERPOLATION':
            diameter = self.diameter
            neuron_flag = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'neuron_flag')
            node_channels = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'node_channels')
            node_length = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'node_length')
            paranodal_length_1 = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                             'paranodal_length_1')
            fiber_type = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'fiber_type')
            self.passive_end_nodes = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,'passive_end_nodes')

            if self.diameter >= 5.643:
                self.delta_z = eval(self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                           'delta_z', 'diameter_greater_or_equal_5.643um'))
            else:
                self.delta_z = eval(self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                           'delta_z', 'diameter_less_5.643um'))
            delta_z = self.delta_z
            paranodal_length_2 = eval(self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode,
                                             'paranodal_length_2'))
            nl = -0.4749 * self.diameter ** 2 + 16.85 * self.diameter - 0.7648
            nodeD = 0.01093 * self.diameter ** 2 + 0.1008 * self.diameter + 1.099
            paraD1 = nodeD
            paraD2 = 0.02361 * self.diameter ** 2 + 0.3673 * self.diameter + 0.7122
            axonD = paraD2

        # Determine number of axonnodes
        if neuron_flag == 2:
            self.axonnodes = int(1 + (n_fiber_coords - 1) / 11)
        elif neuron_flag == 3:
            self.axonnodes = int(n_fiber_coords)
            length = self.delta_z*self.axonnodes

        # Determine starting voltage of fiber type
        if fiber_type == 1:
            self.v_init = -88.3
        elif fiber_type == 2:
            self.v_init = -80
        elif fiber_type == 3:
            channels_type = self.search(Config.FIBER_Z, 'fiber_type_parameters', self.fiber_mode, 'channels_type')
            v_init_c_fibers = [-60, -55, -82, -48] # v_rest for  Sundt, Tigerholm, Rattay and Aberham, and Schild C-Fiber models
            self.v_init = v_init_c_fibers[channels_type - 1]

        # Create fiber sections
        if self.myelination:
            self.createMyelinatedFiber(node_channels, self.axonnodes, self.diameter, self.temperature, axonD,
                                       nodeD, paraD1, paraD2, self.delta_z, paranodal_length_2,
                                       nl, self.passive_end_nodes)
        elif not self.myelination:
            self.createUnmyelinatedFiber(self.diameter, length, c_fiber_model_type=channels_type, celsius=self.temperature,
                                         delta_z=self.delta_z, passive_end_nodes=self.passive_end_nodes)

        return self

    def createMyelinatedFiber(self, node_channels, axonnodes, fiberD, celsius, axonD, nodeD, paraD1, paraD2, deltaz,
                              paralength2, nl, passive_end_nodes):
        def create_MYSA(i, fiberD, paralength1, rhoa, paraD1, e_pas_Vrest, Rpn1, mycm, mygm, nl):
            """
            Create a MYSA segment for MRG_DISCRETE fiber type
            """
            MYSA = Section(name='MYSA ' + str(i))
            MYSA.nseg = 1
            MYSA.diam = fiberD
            MYSA.L = paralength1
            MYSA.Ra = rhoa * (1 / (paraD1 / fiberD) ** 2) / 10000
            MYSA.cm = 2 * paraD1 / fiberD
            MYSA.insert('pas')
            MYSA.g_pas = 0.001 * paraD1 / fiberD
            MYSA.e_pas = e_pas_Vrest

            MYSA.insert('extracellular')
            MYSA.xraxial[0] = Rpn1
            MYSA.xc[0] = mycm / (nl * 2)  # short circuit
            MYSA.xg[0] = mygm / (nl * 2)  # short circuit

            return MYSA

        def create_FLUT(i, fiberD, paralength2, rhoa, paraD2, e_pas_Vrest, Rpn2, mycm, mygm, nl):
            """
            Create a FLUT segment for MRG_DISCRETE fiber type
            """
            FLUT = Section(name='FLUT ' + str(i))
            FLUT.nseg = 1
            FLUT.diam = fiberD
            FLUT.L = paralength2
            FLUT.Ra = rhoa * (1 / (paraD2 / fiberD) ** 2) / 10000
            FLUT.cm = 2 * paraD2 / fiberD
            FLUT.insert('pas')
            FLUT.g_pas = 0.0001 * paraD2 / fiberD
            FLUT.e_pas = e_pas_Vrest

            FLUT.insert('extracellular')
            FLUT.xraxial[0] = Rpn2
            FLUT.xc[0] = mycm / (nl * 2)  # short circuit
            FLUT.xg[0] = mygm / (nl * 2)  # short circuit

            return FLUT

        def create_STIN(i, fiberD, interlength, rhoa, axonD, e_pas_Vrest, Rpx, mycm, mygm, nl):
            """
            Create a STIN segment for MRG_DISCRETE fiber type
            """
            STIN = Section(name='STIN ' + str(i))
            STIN.nseg = 1
            STIN.diam = fiberD
            STIN.L = interlength
            STIN.Ra = rhoa * (1 / (axonD / fiberD) ** 2) / 10000
            STIN.cm = 2 * axonD / fiberD
            STIN.insert('pas')
            STIN.g_pas = 0.0001 * axonD / fiberD
            STIN.e_pas = e_pas_Vrest

            STIN.insert('extracellular')
            STIN.xraxial[0] = Rpx
            STIN.xc[0] = mycm / (nl * 2)  # short circuit
            STIN.xg[0] = mygm / (nl * 2)  # short circuit

            return STIN

        def create_node(index, nodeD, nodelength, rhoa, mycm, mygm, rhoe, passive, axonnodes, node_channels, nl,
                        Rpn0, celsius):
            """
            Create a node of Ranvier for MRG_DISCRETE fiber type
            """
            node = Section(name='node ' + str(index))
            node.nseg = 1
            node.diam = nodeD
            node.L = nodelength
            node.Ra = rhoa / 10000

            if passive and (index == 0 or index == axonnodes - 1):
                node.cm = 2
                node.insert('pas')
                node.g_pas = 0.0001
                node.e_pas = -70
                node.insert('extracellular')
                node.xc[0] = mycm / (nl * 2)  # short circuit
                node.xg[0] = mygm / (nl * 2)  # short circuit

            else:
                if node_channels == 0:
                    node.cm = 2
                    node.insert('axnode_myel')
                elif node_channels == 1:
                    node.cm = 1.149452367  # [uF/cm^2] specific membrane capacitance (Schild 1994, A-type)
                    F = 96500  # [C/mole] Faraday'node Constant from Schild 1994
                    R = 8314  # [J/(kg*mole*K)] Gas Constant from Schild 1994

                    # Based on Schild 1994 ion channels
                    node.insert('leakSchild')
                    node.insert('naf')
                    node.insert('nas')
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
                    node.insert('NakpumpSchild')

                    node.L_caintscale = node.L
                    node.nseg_caintscale = node.nseg
                    node.L_caextscale = node.L
                    node.nseg_caextscale = node.nseg

                    # Ionic concentrations
                    node.cao0_ca_ion = 2.0  # [mM] Initial Cao Concentration
                    node.cai0_ca_ion = 0.000117  # [mM] Initial Cai Concentrations
                    node.ko = 5.4  # [mM] External K Concentration
                    node.ki = 145.0  # [mM] Internal K Concentration
                    kstyle = ion_style("k_ion", 1, 2, 0, 0, 0)  # Allows ek to be calculated manually
                    node.ek = ((R * (celsius + 273.15)) / F) * np.log10(
                        node.ko / node.ki)  # Manual Calculation of ek in order to use Schild F and R values

                    node.nao = 154  # [mM] External Na Concentration
                    node.nai = 8.9  # [mM] Internal Na Concentration
                    nastyle = ion_style("na_ion", 1, 2, 0, 0, 0)  # Allows ena to be calculated manually
                    node.ena = ((R * (celsius + 273.15)) / F) * np.log10(
                        node.nao / node.nai)  # Manual Calculation of ena in order to use Schild F and R values

                    node.gbar_naf = 3  # 0.072503919  # NOTE: Does not conduct with original Schild 1994 value of 0.072503919; does not conduct at 0.5, but does at 1; increased to MRG value of 3
                    node.shiftnaf_naf = 0  # [mV]
                    node.gbar_nas = 3.53678E-07
                    node.shiftnas_nas = 0
                    node.gbar_kd = 0.000194523
                    node.shiftkd_kd = 0
                    node.gbar_ka = 0.001237872
                    node.shiftka_ka = 0
                    node.gbar_kds = 0.000353678
                    node.shiftkds_kds = 0
                    node.gbar_kca = 0.00022989
                    node.gbar_can = 3.53678E-05
                    node.shiftcan_can = 0
                    node.gbar_cat = 1.23787E-05
                    node.shiftcan_cat = 0
                    node.gbna_leak = 1.14945E-05

                node.insert('extracellular')
                node.xraxial[0] = Rpn0
                node.xc[0] = 0  # short circuit
                node.xg[0] = 1e10  # short circuit

            return node

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
            new_node = create_node(i, nodeD, nodelength, rhoa, mycm, mygm, rhoe, passive_end_nodes, axonnodes,
                                    node_channels, nl, Rpn0, celsius)
            self.node.append(new_node)
        for i in range(0, paranodes1):
            new_MYSA = create_MYSA(i, fiberD, paralength1, rhoa, paraD1, e_pas_Vrest, Rpn1, mycm, mygm, nl)
            self.MYSA.append(new_MYSA)
        for i in range(0, paranodes2):
            new_FLUT = create_FLUT(i, fiberD, paralength2, rhoa, paraD2, e_pas_Vrest, Rpn2, mycm, mygm, nl)
            self.FLUT.append(new_FLUT)
        for i in range(0, axoninter):
            new_STIN = create_STIN(i, fiberD, interlength, rhoa, axonD, e_pas_Vrest, Rpx, mycm, mygm, nl)
            self.STIN.append(new_STIN)

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

    def createUnmyelinatedFiber(self, fiberD=6, length=21, c_fiber_model_type=1, celsius=37, delta_z=50/6, insert97na=0,
                     conductances97=0, passive_end_nodes=0):
        """
        Create a list of Neuron Sections for an unmyelinated fiber
        """
        nsegments = int(length/delta_z)

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
                    node.insert('koiTiger')

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

                node.insert('extracellular')
                node.xc[0] = 0  # short circuit
                node.xg[0] = 1e10  # short circuit

        for i in range(0, nsegments-1):
            self.sec[i+1].connect(self.sec[i])

        return self

    def finite_amplitudes(self, stimulation, saving, recording, start_time, amps):
        time_total = 0
        for amp_ind, amp in enumerate(amps):
            print(f'Running amp {amp_ind} of {len(amps)}: {amp} nA')
            self.run(amp, stimulation, recording, saving)
            time_individual = time.time() - start_time - time_total
            saving.saveVariables(self, recording, stimulation.dt, amp_ind)
            saving.saveActivation(self, amp_ind)
            saving.saveRuntime(self, time_individual, amp_ind)

            time_total += time_individual
            recording.reset()

    def findThresh(self, stimulation, saving, recording, find_block_thresh=False):
        # self.run(-1, stimulation, recording, find_block_thresh, saving=saving)
        # return

        bounds_search_mode = self.search(Config.SIM, "protocol", "bounds_search", "mode")
        if bounds_search_mode == 'PERCENT_INCREMENT':  # relative increment (increase bound by a certain percentage of the previous value)
            increment_flag = SearchAmplitudeIncrementMode.PERCENT_INCREMENT.value
            step = self.search(Config.SIM, "protocol", "bounds_search", "step")
            rel_increment = round(step / 100, 4)
        elif bounds_search_mode == 'ABSOLUTE_INCREMENT':  # absolute increment (increase bound by a a certain amount + previous value)
            increment_flag = SearchAmplitudeIncrementMode.ABSOLUTE_INCREMENT.value
            step = self.search(Config.SIM, "protocol", "bounds_search", "step")
            abs_increment = round(step, 4)

        termination_criteria_mode = self.search(Config.SIM, "protocol", "termination_criteria", "mode")
        if termination_criteria_mode == 'ABSOLUTE_DIFFERENCE':
            termination_flag = TerminationCriteriaMode.ABSOLUTE_DIFFERENCE.value
            res = self.search(Config.SIM, "protocol", "termination_criteria", "tolerance")
            abs_thresh_resoln = round(res, 4)
        elif termination_criteria_mode == 'PERCENT_DIFFERENCE':
            termination_flag = TerminationCriteriaMode.PERCENT_DIFFERENCE.value
            res = self.search(Config.SIM, "protocol", "termination_criteria", "percent")
            rel_thresh_resoln = round(res / 100, 4)

        stimamp_top = self.search(Config.SIM, 'protocol', 'bounds_search', 'top')
        stimamp_bottom = self.search(Config.SIM, 'protocol', 'bounds_search', 'bottom')

        check_top_flag = 0  # 0 for upper-bound not yet found, value changes to 1 when the upper-bound is found
        check_bottom_flag = 0  # 0 for lower-bound not yet found, value changes to 1 when the lower-bound is found
        # enter binary search when both are found

        iter = 1
        while True:
            if check_top_flag == 0:
                print("Running stimamp_top = {:.6f}".format(stimamp_top))
                self.run(stimamp_top, stimulation, recording, find_block_thresh)

                if self.n_aps == 0:
                    if find_block_thresh == NeuronRunMode.ACTIVATION_THRESHOLD.value:
                        print(
                            "ERROR: Initial stimamp_top value does not elicit an AP - need to increase its magnitude and/or increase tstop to detect evoked AP")
                    else:
                        print(
                            "WARNING: Initial stimamp_top value does not block - need to increase its magnitude and/or increase tstop to block test pulse evoked AP")
                    if increment_flag == SearchAmplitudeIncrementMode.ABSOLUTE_INCREMENT.value:
                        stimamp_top = stimamp_top + abs_increment
                    elif increment_flag == SearchAmplitudeIncrementMode.PERCENT_INCREMENT.value:
                        stimamp_top = stimamp_top * (1 + rel_increment)
                else:
                    check_top_flag = 1

            if check_bottom_flag == 0:
                print("Running stimamp_bottom = {:.6f}".format(stimamp_bottom))
                self.run(stimamp_bottom, stimulation, recording, find_block_thresh)

                if self.n_aps != 0:
                    if find_block_thresh == NeuronRunMode.ACTIVATION_THRESHOLD.value:
                        print(
                            "ERROR: Initial stimamp_bottom value elicits an AP - need to decrease its magnitude and/or increase tstop to detect block test pulses")
                    else:
                        print(
                            "WARNING: Initial stimamp_bottom value blocks - need to decrease its magnitude and/or increase tstop to detect test pulse evoked AP")
                    if increment_flag == SearchAmplitudeIncrementMode.ABSOLUTE_INCREMENT.value:
                        stimamp_bottom = stimamp_bottom - abs_increment
                    elif increment_flag == SearchAmplitudeIncrementMode.PERCENT_INCREMENT.value:
                        stimamp_bottom = stimamp_bottom * (1 - rel_increment)
                else:
                    check_bottom_flag = 1

            if check_bottom_flag == 1 and check_top_flag == 1:
                print('Bounds set - entering binary search')
                break

            iter += 1

            if iter >= 100:
                print("maximum number of bounds searching steps reached. breaking.")
                quit()

        # enter binary search
        while True:
            stimamp_prev = stimamp_top

            stimamp = (stimamp_bottom + stimamp_top) / 2
            print("stimamp_bottom = {:.6f}      stimamp_top = {:.6f}".format(stimamp_bottom, stimamp_top))
            print("Running stimamp: {:.6f}".format(stimamp))
            self.run(stimamp, stimulation, recording, find_block_thresh)

            if termination_flag == TerminationCriteriaMode.PERCENT_DIFFERENCE.value:
                thresh_resoln = abs(rel_thresh_resoln)
                tolerance = abs((stimamp_bottom - stimamp_top) / stimamp_top)
            elif termination_flag == TerminationCriteriaMode.ABSOLUTE_DIFFERENCE.value:
                thresh_resoln = abs(abs_thresh_resoln)
                tolerance = abs(stimamp_bottom - stimamp_top)

            if tolerance < thresh_resoln:
                if self.n_aps < 1:
                    stimamp = stimamp_prev
                print(
                    "Done searching! stimamp: {:.6f} mA for extracellular and nA for intracellular (check flag_whichstim)\n".format(
                        stimamp))
                self.run(stimamp, stimulation, recording, find_block_thresh, saving=saving)
                saving.saveThresh(self, stimamp)
                break
            elif self.n_aps >= 1:
                stimamp_top = stimamp
            elif self.n_aps < 1:
                stimamp_bottom = stimamp
        return

    def submit(self, stimulation, saving, recording, start_time):
        # determine protocol
        protocol_mode = self.search(Config.SIM, 'protocol', 'mode')
        if protocol_mode != 'FINITE_AMPLITUDES':
            find_thresh = True
            if protocol_mode == 'BLOCK_THRESHOLD':
                find_block_thresh = NeuronRunMode.BLOCK_THRESHOLD.value
            elif protocol_mode == 'ACTIVATION_THRESHOLD':
                find_block_thresh = NeuronRunMode.ACTIVATION_THRESHOLD.value
        elif protocol_mode == 'FINITE_AMPLITUDES':
            find_thresh = False
            find_block_thresh = False
            amps = self.search(Config.SIM, 'protocol', 'amplitudes')

        if find_thresh:
            self.findThresh(stimulation, saving, recording, find_block_thresh)
            saving.saveVariables(self, recording, stimulation.dt)
            time_individual = time.time()-start_time
            saving.saveRuntime(self, time_individual)

        else:
            self.finite_amplitudes(stimulation, saving, recording, start_time, amps)

    def run(self, stimamp, stimulation, recording, find_block_thresh=False, saving=None):
        """
        Run a simulation for a single stimulation amplitude
        """
        def balance(fiber):
            # Balance membrane currents for Tigerholm model
            Vrest = -55
            for s in fiber.sec:
                if (-(s.ina_nattxs + s.ina_nav1p9 + s.ina_nav1p8 + s.ina_h + s.ina_nakpump) / (Vrest - s.ena)) < 0:
                    s.pumpina_extrapump = -(s.ina_nattxs + s.ina_nav1p9 + s.ina_nav1p8 + s.ina_h + s.ina_nakpump)
                else:
                    s.gnaleak_leak = -(s.ina_nattxs + s.ina_nav1p9 + s.ina_nav1p8 + s.ina_h + s.ina_nakpump) / (
                                Vrest - s.ena)

                if (-(s.ik_ks + s.ik_kf + s.ik_h + s.ik_kdrTiger + s.ik_nakpump + s.ik_kna) / (Vrest - s.ek)) < 0:
                    s.pumpik_extrapump = -(s.ik_ks + s.ik_kf + s.ik_h + s.ik_kdrTiger + s.ik_nakpump + s.ik_kna)
                else:
                    s.gkleak_leak = -(s.ik_ks + s.ik_kf + s.ik_h + s.ik_kdrTiger + s.ik_nakpump + s.ik_kna) / (
                                Vrest - s.ek)

        def steady_state(fiber, stim_dt):
            # Allow system to reach steady-state by using a large dt before simulation
            t_initSS = fiber.search(Config.SIM, 'protocol', 'initSS')
            dt_initSS = fiber.search(Config.SIM, 'protocol', 'dt_initSS')
            h.t = t_initSS      # Start before t=0
            h.dt = dt_initSS    # Large dt
            while (h.t <= -dt_initSS):
                h.fadvance()
            h.dt = stim_dt
            h.t = 0
            h.fcurrent()
            h.frecord_init()

        # If saving variables, record variables
        if saving is not None:
            if saving.time_vm or saving.space_vm:
                recording.record_vm(self)
            if saving.time_gating or saving.space_gating:
                recording.record_gating(self)
            if saving.istim:
                recording.record_istim(stimulation.istim)
            if saving.ap_end_times:
                recording.record_ap_end_times(self, saving.ap_end_inds, saving.ap_end_thresh)

        h.finitialize(self.v_init)          # Initialize the simulation
        if self.fiber_mode == 'TIGERHOLM':  # Balance membrane currents if Tigerholm
            balance(self)

        stimulation.initialize_extracellular(self)  # Set extracellular stimulation at each segment to zero
        steady_state(self, stimulation.dt)          # Allow system to reach steady-state before simulation
        h.celsius = self.temperature                # Set simulation temperature

        # Set up APcount
        recording.record_ap(self)

        ############    Begin simulation     ############
        n_tsteps = len(stimulation.waveform)
        for i in range(0, n_tsteps):
            if h.t > stimulation.tstop:
                break
            amp = stimulation.waveform[i]
            scaled_stim = [stimamp*amp*x for x in stimulation.potentials]
            stimulation.update_extracellular(self, scaled_stim)

            h.fadvance()
        ############    Done with simulation     ############

        # Insert vectors of 0's for gating parameters at passive end nodes
        if saving is not None:
            if saving.time_gating or saving.space_gating:
                recording.record_gating(self, fix_passive=True)

        # Check if APs occurred
        self.n_aps = recording.ap_checker(self, find_block_thresh)

        if saving is None:
            print(f'{int(self.n_aps)} AP(s) detected')
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
