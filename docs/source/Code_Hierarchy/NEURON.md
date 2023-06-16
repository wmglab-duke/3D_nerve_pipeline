# NEURON Files

## NEURON simulations
All NEURON simulations are handled using the `wmglab-neuron` PyPI package (Github: _link_).

### Create fiber model

Based on the "fiber_model" (associated by the "fiber/modes" type parameter
 in `<sim_index>.json` ([Sim Parameters](../JSON/JSON_parameters/sim))),
`run_controls.py` builds an instance of the fiber model using the `wmglab-neuron` package.
For all fiber types, the segments created and connected in NEURON have lengths that correspond to
the coordinates of the input potentials.

### Intracellular stimulus

For simulations of block threshold, an intracellular test pulse is
delivered at one end of the fiber to test if the cuff electrode (i.e.,
placed between the intracellular stimulus and the site of detecting
action potentials) is blocking action potentials ([Simulation Protocols](../Running_ASCENT/Info.md#simulation-protocols)). The intracellular
stimulation parameters are defined in **_Sim_** and are defined as
parameters in NEURON within the `run_controls.py` file. The parameters in
**_Sim_** control the pulse delay, pulse width, pulse repetition
frequency, pulse amplitude, and node/section index of the intracellular
stimulus ([Sim Parameters](../JSON/JSON_parameters/sim)). For simulating activation thresholds, the intracellular
stimulation amplitude should be set to zero.

### Extracellular stimulus

To simulate response of individual fibers to electrical stimulation, we
use NEURON’s extracellular mechanisms in the `wmglab-neuron` package
to apply the electric potential from COMSOL at each segment of the
cable model as a time-varying signal.
`run_controls.py` loads in the stimulation waveform from a `n_sim’s` `data/inputs/`.
The saved stimulation waveform is unscaled,
meaning the maximum current magnitude at any timestep is +/-1.
Analogously, `run_control.py` reads in the potentials for the fiber being simulated
from `data/inputs/`.

### Recording

The NEURON simulation code contains functionality ready to record and
save to file the values of state variables at discrete spatial locations
for all times and/or at discrete times for all spatial locations (i.e.,
nodes of Ranvier for myelinated fibers or sections for unmyelinated
fibers) for applied extracellular potential, intracellular stimulation
amplitude, transmembrane potential, and gating parameters using
the `wmglab-neuron` package. The recording tools are particularly useful for
generating data to troubleshoot and visualize simulations.

### Running a simulation & different protocols

The `wmglab-neuron` package is responsible for simulating the response of the
model fiber to intracellular and extracellular stimulation. Similarly, the package handles
various simulation protocols ([Simulation Protocols](../Running_ASCENT/Info.md#simulation-protocols)).

### Save outputs to file

At the end of the NEURON simulation, the program saves state variables
as indicated with saveflags, CPU time, and threshold values. Output
files are saved to the `data/outputs/` directory within its `n_sim` folder.

## run_controls.py

The `run_controls.py` file controls the parameters and simulation protocol for
modeling fiber response to electrical stimulation in NEURON and is
automatically populated based on parameters in **_Model_** and
**_Sim_**.
Parameters defined in `<sim_index>.py` and `<model_index>.json` span the categories of: environment
(i.e., temperature from **_Model_**), simulation time (i.e., time step,
duration of simulation from **_Sim_**), fiber parameters (i.e., flags
for fiber geometry and channels, number of fiber nodes from **_Model_**,
**_Sim_**, and `config/system/fiber_z.json`), intracellular stimulation
(i.e., delay from start of simulation, amplitude, pulse duration, pulse
repetition frequency from **_Sim_**), extracellular stimulation (i.e.,
path to waveform file in `n_sim/` folder which is always
`data/inputs/waveform.dat`), flags to define the model parameters that
should be recorded (i.e., Vm(t), Gating(t), Vm(x), Gating(x) from
**_Sim_**), the locations at which to record the parameters (nodes of
Ranvier for myelinated axons from **_Sim_**), and parameters for the
bisection search for thresholds (i.e., activation or block protocol,
initial upper and lower bounds on the stimulation amplitude for the
bisection search, and threshold resolution for the bisection search from
**_Sim_**). The `run_controls.py` file takes inputs of the **_Sim_** directory, `n_sim/` directory,
and an exception configuration. The file automatically loads the `fiber_z.json` configuration
file which contains all associated flags, parameters, and rules for
defining a fiber’s geometry and channel mechanisms in NEURON, all of which is handled in the
`wmglab-neuron` package.
