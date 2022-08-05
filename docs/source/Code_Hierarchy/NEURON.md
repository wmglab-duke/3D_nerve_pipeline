# NEURON Files

## Overview

All of ASCENT's NEURON processes have been refactored from HOC into Python.
There are five new Python files that work together to create and simulate electrical stimulation of a fiber: `fiber.py`, `run_controls.py`, `recording.py`, `stimulation.py`, and `saving.py`.
During the pipeline, each n_sim gets an instance of the Fiber class in `Fiber.py`, which is saved to `submit/n_sims/<n_sim_index>/fiber.obj`.
During NEURON job submissions, this Fiber object is loaded and builds its own NEURON sections according to the given `<n_sim_index>.json` file, runs its own binary search for thresholds (i.e. activation or block protocol) or amplitudes (i.e. finite amplitudes protocol), and runs individual run simulations.
The rest of the classes are used to coordinate all program operations (`run_controls.py`), manage stimulation (`stimulation.py`), record run simulation data (`recording.py`), and output simulation data to file (`stimulation.py`).

## Create fiber model

In `run_controls.py`, the simulation Fiber object is loaded from the `n_sim/<n_sim_index>` directory.
The `Fiber.generate()` method builds the appropriate NEURON Sections based on the user-specified fiber model type. `Fiber.generate()` loads the `fiber_z.json` configuration for the given fiber type containing all associated flags, parameters, and rules for
defining a fiber’s geometry and channel mechanisms in NEURON. `Fiber.generate()` determines the number of axon nodes (nodes of Ranvier for myelinated fibers and axon segments for unmyelinated fibers) and resting membrane voltage before calling `Fiber.createMyelinatedFiber()` or `Fiber.createUnmyelinatedFiber()` to create the actual fiber sections.
For all fiber types, the segments created and connected in NEURON have lengths that correspond to the coordinates of the input potentials.
If the user-specified fiber model type is myelinated, `Fiber.createMyelinatedFiber()` creates a list of NEURON sections for each of the four axon segment types (node of Ranvier, FLUT, MYSA, and STIN) and connects the Sections.
`Fiber.createMyelinatedFiber()` contains submethods `create_node()`, `create_MYSA()`, `create_FLUT()`, and `create_STIN()`, which returns individual NEURON sections for node of Ranvier, MYSA, FLUT, and STIN axon segments, respectively.
If the user specifies a fiber model type that is unmyelinated, `Fiber.createUnmyelinatedFiber()` creates a single list of NEURON sections for unmyelinated axon segment types and connects each segment sequentially.

## Intracellular stimulus

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

## Extracellular stimulus

To simulate response of individual fibers to electrical stimulation, we
use NEURON’s extracellular mechanisms to apply the electric potential
from COMSOL at each segment of the cable model as a time-varying signal.
We load in the stimulation waveform, as well as simulation time variables (i.e., time step, duration of simulation), from a `n_sim’s` `data/inputs/`
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

## Recording

The NEURON simulation code contains functionality ready to record and
save to file the values of state variables at discrete spatial locations
for all times and/or at discrete times for all spatial locations (i.e.,
nodes of Ranvier for myelinated fibers or sections for unmyelinated
fibers) for applied extracellular potential, intracellular stimulation
amplitude, transmembrane potential, and gating parameters using
the Recording class in `recording.py`. The recording tools are particularly useful for
generating data to troubleshoot and visualize simulations.

## Protocols

The method `Fiber.submit()` in `Fiber.py` determines protocol from **_Sim_** and makes calls to the appropriate functions.
The method `Fiber.findThresh()` in `Fiber.py` performs a binary search for activation and
block thresholds. The method `Fiber.finite_amplitudes()` in `Fiber.py` iterates over each of the amplitudes specified by the user in
**_Sim_**. Both of these methods make calls to `Fiber.run()` in `Fiber.py` which simulates an individual run simulation for a given stimulation amplitude.

See documentation on Simulation Protocols for more information ([Simulation Protocols](../Running_ASCENT/Info.md#simulation-protocols)).

## Running an individual run simulation

The method `Fiber.run()` in `Fiber.py` is responsible for simulating the response of the
model fiber to intracellular and extracellular stimulation. Before the
simulation starts, the procedure adds action potential counters to look
for a rise above a threshold transmembrane potential.

So that each fiber reaches a steady-state before the simulation starts,
the `Fiber.run()` method contains the submethod `steady_state()` that initializes the fiber by stepping through large
time steps with no extracellular potential applied to each compartment.
`Fiber.run()` then loops over each time step, and, while updating the value of
extracellular potential at each fiber segment, records the values of
flagged state variables as necessary.

After `Fiber.run()` loops over all time steps, if the user is
searching for threshold current amplitudes, the method evaluates if the
extracellular stimulation amplitude was above or below threshold, as
indicated by the presence or absence of an action potential for
activation and block thresholds, respectively.

## Save outputs to file

At the end of the NEURON simulation, the program saves state variables
as indicated in **_Sim_**, CPU time, and threshold values (activation and block protocols) or number of action potentials (finite amplitudes protocol). Output
files are saved to the `data/outputs/` directory within its `n_sim` folder.

## fiber.py

The `fiber.py` file contains the Fiber class, which is at the heart of all of ASCENT's NEURON processes. An instance of the Fiber class is first created in `Simulation.n_sim_setup()` during the ASCENT Pipeline.
At this point in the pipeline, an instance of the Fiber class is created, configured with the specific n_sim's **_Model_**, **_Sim_**, and `fiber_z.json`. `Simulation.n_sim_setup()` calls `Fiber.inherit()`, which initializes the object with universal attributes for all fibers in a given n_sim (i.e., fiber diameter, min, max, offset, seed from `fibers/z_parameters` in **_Sim_**, fiber model type from **_Sim_**, and model temperature from **_Model_**).
`Simulation.n_sim_setup()` then saves this instance of the Fiber class to `n_sim/<n_sim_index>/fiber.obj` for the appropriate directory.
During NEURON job submissions, this Fiber object is loaded from `n_sim/<n_sim_index>/fiber.obj`. The Fiber class reads from **_Sim_** and `fiber_z.json` to determine all associated flags, parameters, and rules for
defining a fiber’s geometry and channel mechanisms in NEURON. The `Fiber.generate()` method builds the appropriate NEURON Sections based on the user-specified fiber model type. `Fiber.generate()` loads the `fiber_z.json` configuration for the given fiber type containing all associated flags, parameters, and rules for
defining a fiber’s geometry and channel mechanisms in NEURON. `Fiber.generate()` determines the number of axon nodes (nodes of Ranvier for myelinated fibers and axon segments for unmyelinated fibers) and resting membrane voltage before calling `Fiber.createMyelinatedFiber()` or `Fiber.createUnmyelinatedFiber()` to create the actual fiber sections.
For all fiber types, the segments created and connected in NEURON have lengths that correspond to the coordinates of the input potentials.
If the user-specified fiber model type is myelinated, `Fiber.createMyelinatedFiber()` creates a list of NEURON sections for each of the four axon segment types (node of Ranvier, FLUT, MYSA, and STIN) and connects the Sections.
`Fiber.createMyelinatedFiber()` contains submethods `create_node()`, `create_MYSA()`, `create_FLUT()`, and `create_STIN()`, which returns individual NEURON sections for node of Ranvier, MYSA, FLUT, and STIN axon segments, respectively.
If the user specifies a fiber model type that is unmyelinated, `Fiber.createUnmyelinatedFiber()` creates a single list of NEURON sections for unmyelinated axon segment types and connects each segment sequentially.
The method `Fiber.submit()` in `Fiber.py` determines protocol from **_Sim_** and makes calls to the appropriate functions.
The method `Fiber.findThresh()` in `Fiber.py` performs a binary search for activation and
block thresholds. The method `Fiber.finite_amplitudes()` in `Fiber.py` iterates over each of the amplitudes specified by the user in
**_Sim_**. Both of these methods make calls to `Fiber.run()` in `Fiber.py` which simulates an individual run simulation for a given stimulation amplitude.

## run_controls.py

The `run_controls.py` file coordinates all program operations to create a
biophysically realistic discrete cable fiber model, simulate the fiber’s
response to extracellular and intracellular stimulation, and record the
response of the fiber. For each fiber simulation, an instance of the Fiber
object (`src/core/fiber.py`) is loaded from the `n_sim/<n_sim_index>/` directory.
Instances of the Recording class (`recording.py`), Stimulation class (`stimulation.py`), and Saving class (`saving.py`)
are created and initialized from parameters in **_Sim_**.

## recording.py

The Recording class in `recording.py` is responsible for recording all data during a single run simulation.
The Recording class uses Python lists, NEURON Vector() objects, and NEURON APCount() objects to record spacial distribution, time, and user-specified state variables.
For all run simulations, the Recording class records action potentials (`Recording.record_ap()`) at all segments (unmyelinated) or nodes of Ranvier (myelinated) and
checks at the end of the simulation for the number of action potentials that occurred at the location specified by the user in **_Sim_** (`Recording.ap_checker()`).
Depending on the state variables the user chooses to save in **_Sim_**, for a given run simulation, the Recording class records the following:
membrane potential (`Recording.record_vm()`), channel gating parameters (`Recording.record_gating()`), applied intracellular stimulation (`Recording.record_istim()`), and further information about action potentials (`Recording.record_ap_end_times()`).
The Recording class can also reset its recording data for subsequent run simulations (`Recording.reset()`), which is required for the finite amplitude protocol, where multiple runs require recorded data to be saved to file.
For activation and block thresholds, these variables are only saved during a final run simulation at the threshold current amplitude. For finite amplitude protocols, these variables are recorded for every run simulation.

## saving.py

The `saving.py` file contains the Saving class. The Saving class is responsible for outputting all simulation data to file.
The method `Saving.inherit()` sets the Saving class attributes according to the parameters in **_Sim_**, indicating what state variables to output to file and at what times and locations to save said state variables, as well as the path to the output directory for the given n_sim.
Outputs are saved to `<n_sim_index>/data/outputs/` as text files.
Depending on **_Sim_**, `Saving.saveVariables()` outputs state variables at discrete times and/or locations to file. `Saving.saveVariables()` contains a submethod `create_header()`, which creates a list of strings to be used as the first line of the outputted text files.
Analogously, `Saving.saveRuntime()` outputs the CPU time for an individual simulation to file, if it is indicated by the user in **_Sim_**.
For activation and block threshold protocols, the method `Saving.saveThresh()` outputs threshold current amplitudes to file.
For the finite amplitude protocol, the method `Saving.saveActivation()` outputs the number of action potentials that occurred at the location specified in **_Sim_** to file.

## stimulation.py

The `stimulation.py` file contains the Stimulation class. The Stimulation class is responsible for all data related to extracellular stimulation (fiber potentials and waveform) and intracellular stimulation.
At the beginning of all NEURON processes, `Stimlation.load_potentials()` loads fiber potentials and `Stimulation.load_waveform()` loads the waveform from file (`inner#_fiber#.dat` and `waveform.dat`, respectively, from `<n_sim_index>/data/inputs/`).
`Stimulation.load_waveform()` also determines simulation time (i.e., time step,
duration of simulation) from the `waveform.dat` file as well.
Similarly, `Stimulation.apply_intracellular()` creates an instance of the NEURON trainIClamp object.
Before each simulation in `Fiber.run()`, `Stimulation.initialize_extracellular()` initializes all extracellular stimulation,
and at each time step of the simulation, `Stimulation.update_extracellular()` updates the applied extracellular stimulation at each section of the fiber.
