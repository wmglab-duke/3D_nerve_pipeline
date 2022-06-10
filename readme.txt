How to reproduce threshold values:

- set up env.json (might not be necessary)
- make sure run config exists in ascent/config/user/runs/
- make sure potentials from COMSOL are in the path: ascent/samples/<sample_index>/models/<model_index>/sims/<sim_index>/n_sims/<nsim_index>/data/inputs/
- make sure no sim.obj exists for in ascent/samples/<sample_index
>/models/<model_index>/sims/<sim_index>/
- make sure hardcoded paths are correct in run and generate methods in fiber.py (lines 163, 442, 443)
- run "python run pipeline <run_index>" (run_index = 0 for beta task)
