import pandas as pd
import numpy as np
import pickle
from neuron import h
import bokeh.plotting as plt

h.load_file('stdrun.hoc')

Section = h.Section

fiber = pickle.load(open('fiber.obj', 'rb'))
fiber.generate(133)

intracellular_stim = h.trainIClamp(fiber.node[1](0.5))
intracellular_stim.delay = 10
intracellular_stim.PW = 0.1
intracellular_stim.train = 1
intracellular_stim.freq = 100
intracellular_stim.amp = 10

plot_node0 = plt.figure(title='node[0]', x_axis_label='t (ms)', y_axis_label='v (mV)')
v_node0 = h.Vector().record(fiber.node[1](0.5)._ref_v)
t = h.Vector().record(h._ref_t)
ap_times = h.Vector()
apc = h.APCount(fiber.node[1](0.5))
apc.thresh = -30
apc.record(ap_times)

h.finitialize(-80)
h.celsius = 37
h.continuerun(100)

plot_node0.line(t, list(v_node0), line_width=2)
plt.show(plot_node0)

print(list(ap_times))