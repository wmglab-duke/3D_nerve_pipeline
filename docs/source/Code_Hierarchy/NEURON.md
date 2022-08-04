# NEURON Files

OVERVIEW:
All of ASCENT's NEURON processes have been refactored from HOC into Python.
There are five new files that work together to create and simulate electrical stimulation of a fiber: `fiber.py`, `run_controls.py`, `recording.py`, `stimulation.py`, and `saving.py`.

## run_controls.py

The `run_controls.py` file coordinates all program operations to create a
biophysically realistic discrete cable fiber model, simulate the fiber’s
response to extracellular and intracellular stimulation, and record the
response of the fiber. For each fiber simulation, an instance of the Fiber
object (`src/core/fiber.py`) is loaded. Instances of the Recording class (`recording.py`), Stimulation class (`stimulation.py`), and Saving class (`saving.py`)
are created and initialized with the data for the specific fiber simulation.

## recording.py

The NEURON simulation code contains functionality ready to record and
save to file the values of state variables at discrete spatial locations
for all times and/or at discrete times for all spatial locations (i.e.,
nodes of Ranvier for myelinated fibers or sections for unmyelinated
fibers) for applied extracellular potential, intracellular stimulation
amplitude, transmembrane potential, and gating parameters using
the Recording class in `recording.py`. The recording tools are particularly useful for
generating data to troubleshoot and visualize simulations. The Recording class is responsible for recording all data during a single run simulation.
The Recording class uses Python lists, NEURON Vector() objects, and APCount() objects to record spacial distribution, time, and user-specified state variables.
For all runs, the Recording class records action potentials at all segments (unmyelinated) or nodes (myelinated) and
checks at the end of the simulation for the number of action potentials that occurred at user-specified location along the axon.
Depending on the state variables the user chooses to save (**_Sim_**), for a given run simulation, the Recording class records the following:
membrane potential, channel gating parameters, applied intracellular stimulation, CPU time, and further information about action potentials (see **_Sim_**).

## saving.py

The `saving.py` file contains the Saving class. The Saving class is responsible for saving all simulation data to file.
The attributes of the Saving class let ASCENT know what state variables to save, at what times and locations to save said state variables, and where to save the outputs to.
For each fiber simulated in NEURON, outputs are saved to `<n_sim_index>/data/outputs/` as text files.
For protocols `BLOCK_THRESHOLD` and `ACTIVATION_THRESHOLD` (SEE PROTOCOLS), data outputs include threshold
current amplitudes. For `FINITE_AMPLITUDES` protocols, data outputs include the number of action potentials that occurred at a user-specified location along the fiber.
Depending on **_Sim_**, data outputs can also include state variables at discrete times and/or locations.

## stimulation.py

The `stimulation.py` file contains the Stimulation class. The Stimulation class is responsible for all data related to extracellular stimulation (fiber potentials and waveform) and intracellular stimulation.
At the beginning of all NEURON processes, fiber potentials and waveform are loaded in from file (`<n_sim_index>/data/inputs/`), and an instance of the NEURON trainIClamp object is created. Certain waveform parameters (time step and time stop) are read as well.
Before each simulation, the Stimulation class initializes all extracellular stimulation, and at each time step of the simulation, the Stimulation class updates the applied extracellular stimulation all along the fiber length.

### Create fiber model

In `run_controls.py`, the simulation Fiber object builds the appropriate NEURON Sections based on the user-specified fiber model type.
For all fiber types, the segments created and connected in NEURON have lengths that correspond to the coordinates of the input potentials.
If the user-specified fiber model type is myelinated, a list of NEURON sections each of the four axon segment types (Node of Ranvier, FLUT, MYSA, and STIN) are created, and each of those Sections are connected according to their order in literature (FIX THIS FOR SURE).
If the user specifies a fiber model type that is unmyelinated, a single list of NEURON sections is created for unmyelinated axon segment types.
For each NEURON Section, channel mechanisms are inserted from MOD Files (LINK TO MOD FILES).

### Intracellular stimulus

For simulations of block threshold, an intracellular test pulse is
delivered at one end of the fiber to test if the cuff electrode (i.e.,
placed between the intracellular stimulus and the site of detecting
action potentials) is blocking action potentials ([Simulation Protocols](../Running_ASCENT/Info.md#simulation-protocols)). The intracellular
stimulation parameters are defined in **_Sim_**. Intracellular stimulation is applied using the `Stimulation.apply_intracellular()`
method in `stimulation.py` using the parameters in **_Sim_**.
The parameters in **_Sim_** control the pulse delay, pulse width, pulse repetition
frequency, pulse amplitude, and node/section index of the intracellular
stimulus ([Sim Parameters](../JSON/JSON_parameters/sim)). For simulating activation thresholds, the intracellular
stimulation amplitude should be set to zero.

### Extracellular stimulus

To simulate response of individual fibers to electrical stimulation, we
use NEURON’s extracellular mechanisms to apply the electric potential
from COMSOL at each segment of the cable model as a time-varying signal.
We load in the stimulation waveform, as well as the simulation time step and stop time, from a `n_sim’s` `data/inputs/`
directory using the `Stimulation.load_waveform()` method in
`stimulation.py`. The saved stimulation waveform is unscaled,
meaning the maximum current magnitude at any timestep is +/-1.
Analogously, we read in the potentials for the fiber being simulated
from `data/inputs/` using the `Stimulation.load_potentials()` method in
`stimulation.py`. Before each simulation, extracellular stimulation is initialized
using the `Stimulation.initialize_extracellular()` method, and at each time step of
the simulation, the extracellular stimulation is updated at all sections of the fiber
using the `Stimulation.update_extracellular()` method, both of which are methods in
`stimulation.py`.

### Recording

The NEURON simulation code contains functionality ready to record and
save to file the values of state variables at discrete spatial locations
for all times and/or at discrete times for all spatial locations (i.e.,
nodes of Ranvier for myelinated fibers or sections for unmyelinated
fibers) for applied extracellular potential, intracellular stimulation
amplitude, transmembrane potential, and gating parameters using
`Recording.hoc`. The recording tools are particularly useful for
generating data to troubleshoot and visualize simulations.

### RunSim

The method `Fiber.run()` in `Fiber.py` is responsible for simulating the response of the
model fiber to intracellular and extracellular stimulation. Before the
simulation starts, the procedure adds action potential counters to look
for a rise above a threshold transmembrane potential.

So that each fiber reaches a steady-state before the simulation starts,
the `Fiber.run()` method contains a `steady_state()` submethod that initializes the fiber by stepping through large
time steps with no extracellular potential applied to each compartment.
`Fiber.run()` then loops over each time step, and, while updating the value of
extracellular potential at each fiber segment, records the values of
flagged state variables as necessary.

After `Fiber.run()` loops over all time steps, if the user is
searching for threshold current amplitudes, the method evaluates if the
extracellular stimulation amplitude was above or below threshold, as
indicated by the presence or absence of an action potential for
activation and block thresholds, respectively.

### FindThresh

The method `Fiber.findThresh()` in `Fiber.py` performs a binary search for activation and
block thresholds ([Simulation Protocols](../Running_ASCENT/Info.md#simulation-protocols)).

### Save outputs to file

At the end of the NEURON simulation, the program saves state variables
as indicated with saveflags, CPU time, and threshold values. Output
files are saved to the `data/outputs/` directory within its `n_sim` folder.

## NEURON launch.hoc

The `launch.hoc` file defines the parameters and simulation protocol for
modeling fiber response to electrical stimulation in NEURON and is
automatically populated based on parameters in **_Model_** and
**_Sim_**. The `launch.hoc` file is created by the `HocWriter` class.
Parameters defined in `launch.hoc` span the categories of: environment
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
binary search for thresholds (i.e., activation or block protocol,
initial upper and lower bounds on the stimulation amplitude for the
binary search, and threshold resolution for the binary search from
**_Sim_**). The `launch.hoc` file loads `Wrapper.hoc` which calls all NEURON
procedures. The `launch.hoc` file is created by the Python `HocWriter`
class, which takes inputs of the **_Sim_** directory, `n_sim/` directory,
and an exception configuration. When the `HocWriter` class is
instantiated, it automatically loads the `fiber_z.json` configuration
file which contains all associated flags, parameters, and rules for
defining a fiber’s geometry and channel mechanisms in NEURON.
