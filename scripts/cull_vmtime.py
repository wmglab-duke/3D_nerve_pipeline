"""Created on Wed Mar  2 17:26:43 2022.

@author: dpm42
"""
import json
import os

import numpy as np


def get_activation(
    filepath,
    delta_V: float = 60,
    rounding_precision: int = 5,
    plot: bool = False,
    plot_nodes_on_find: bool = False,
    plot_compiled: bool = False,
    absolute_voltage: bool = True,
    n_sim_label_override: str = None,
    save: bool = False,
    subplots=False,
    nodes_only=False,
    sample_override=None,
    delete_vmtime=False,
):
    vm_t_data = np.loadtxt(filepath, skiprows=1)

    V_o = np.mean(vm_t_data[0, 1:])

    # find dt by rounding first timestep
    dt = round(vm_t_data[1, 0] - vm_t_data[0, 0], rounding_precision)

    # initialize value AP time, node (locations), voltages at time
    time, node, voltages = None, None, None

    # loop through and enumerate each timestep
    rows = vm_t_data[:, 1:]
    index = int(len(rows) / 2)
    for i, row in enumerate(rows):
        # get list of node indices that satisfy deflection condition
        found_nodes = np.where(row >= V_o + delta_V)[0]
        # that list contains any elements, set time and node (location), then break out of loop
        if len(found_nodes) > 0:
            time = round(i * dt, rounding_precision)
            node = found_nodes[0]
            voltages = row
            index = i
            break
    outdata = {'time': time, 'node': int(node + 1), 'fiber_node_count': len(vm_t_data[0, 1:])}
    outfile = filepath.replace('.dat', '.json')
    with open(outfile, 'w') as f:
        json.dump(outdata, f)
    os.remove(filepath)


path = r'D:\ASCENT\ascent\out\n_sims'
for nsim in os.listdir(path):
    print(nsim)
    if os.path.isdir(os.path.join(path, nsim, 'data', 'outputs')):
        thispath = os.path.join(path, nsim, 'data', 'outputs')
        for file in os.listdir(thispath):
            if file.endswith('_amp0.dat') and file.startswith('Vm_time_'):
                filepath = os.path.join(path, nsim, 'data', 'outputs', file)
                get_activation(filepath)
