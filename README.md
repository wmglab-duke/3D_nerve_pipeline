<!-- TODO
- Add dependencies and versions of code / libraries / software / etc.
- Add note about config README
 -->

# 3D modeling of peripheral nerve stimulation

# Installation

1. Clone the repository
2. Create a subdirectory called "sipsource" and place the input .sip files there.
3. `pip install -r requirements.txt`
4. Copy `config/system/example_env.json` and rename to `config/system/env.json` edit with the following paths:
    - ASCENT_COMSOL_PATH: Path to COMSOL installation ending in `/multiphysics`
    - ASCENT_JDK_PATH: Path to jdk installation ending in `/bin` (Must be either JDK8 or 11, 11 is recommended)
    - ASCENT_PROJECT_PATH: Root of the repository
    - ASCENT_NSIM_EXPORT_PATH: Not used but must be present
    - SIMPLEWARE_SCANIP_PATH: Path to simpleware installation which contains the executable "console scan ip"
5. Run the installation command `python run install`, with the `--no-conda` flag  if you already have an ASCENT environment (but them must manually `pip install -r requirements.txt`)

# Running

1. Create a configuration file (.e.g., `2LDS5.json`) from `example.json` (place in same directory)
2. Change the "sourcefile" to the correct file (name only, not full path).
3. Change the "project_path" to the directory where you want the files to be created
4. Navigate to the "Scripts" dir
5. `python 3D_pipeline_ds.py 2LDS5.json` (replace with the configuration file you created)
6. To generate thresholds, take the directories from the 7_ascent folder after running and place them in the relevant ASCENT directories (under the "3d" branch). Add `"post_java_only":true` to your ASCENT run config.
Note: Simpleware will only run if there is a display, so on the cluster this must be run in a DCC desktop session.
Note2: If using Imthera cuff, your 3D config must have `"use_nastran": true` in your mesh block. If using Livanova, this must be absent or false. Other cuffs are not tested.
Note3: examples/3D has a mock nerve which can be used for local runs
Note4: When generating a sim json config file for the ascent run, keep in mind that the 3D pipeline is currently hard coded to save the super sampled bases (`ss_bases/`) into the `3D/models/0/sims/3/` folder. Therefore, the sim file you create should be set up to use the super sampled bases from `"source_sim":3`.

# Deformation

1. Use a sample where deform is set to True (Note, currently only supports LivaNova cuff and deform_ratio of 1)
2. Pipeline will exit prior to actually running deformation simulation. You need to check the following:
     - Look through predeform.sip to see that no fascicle jump around. If they do, you must manually correct them.
     - Manually run the deform algorithm. You may need to adjust the contact penalty parameters to get a successful simulation. Even then, you may not succeed. However, if you simulation gets most of the way there, it may still be usable since the ramping is exponential decreasing.
     - Export as postdeform.stl
     - Import into Simpleware predeform file and scroll through to check that everything looks ok.
     - Manually use the code in postdeform.py to generate contour coordinates for the deformed fascicles.
     - Run the pipeline again, it should skip deformation now that the proper files are in place.
     - If there are any topological mismatches a window will pop up and ask you to indicate fascicles which watershed segmentation thinks are separate but are not.
     - Once this completes, the pipeline will run as normal.

Note: Sometimes the mesh from comsol will be non manifold. In these cases, you may need to "fix" the import geometry before exporting contours. However, this will result in a LOT of windows popping up during the watershed. To avoid this, new code should be written that loops through and export all possible contours where the mesh is manifold, skipping the others, then fixes the stl, and exports the remaining contours.

Other note: sometimes you may want to run a new cuff on a deformed nerve without redoing the whole deformation process. In this case, follow all the steps to make a new config as normal, then:
1. Copy the "slides" folder to the relevant new directory (in datanew, e.g., from 2Ldeform to 2Ldeformnewcuff)
2. Make sure 3D pipeline has the following settings: slidegen = True, skipsave = True, rundeform = False
3. Run as normal, the pipeline should use the deformed slides and put the new cuff on the geometry
