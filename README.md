<!-- TODO
- Add dependencies and versions of code / libraries / software / etc.
- Add note about config README
 -->

[![DOI](https://zenodo.org/badge/1149201201.svg)](https://doi.org/10.5281/zenodo.18475037)

# 3D modeling of peripheral nerve stimulation

Note that this pipeline requires basic familiarity with the ASCENT pipeline (https://github.com/wmglab-duke/ascent). Therefore, it is recommended to complete the ASCENT tutorial before using this software.

# Installation

1. Clone the repository
3. `pip install -r requirements.txt`
4. Copy `config/system/example_env.json` and rename to `config/system/env.json` edit with the following paths:
    - ASCENT_COMSOL_PATH: Path to COMSOL installation ending in `/multiphysics`
    - ASCENT_JDK_PATH: Path to jdk installation ending in `/bin` (Must be either JDK8 or 11, 11 is recommended)
    - ASCENT_PROJECT_PATH: Root of the repository
    - ASCENT_NSIM_EXPORT_PATH: Not used but must be present
    - SIMPLEWARE_SCANIP_PATH: Path to simpleware installation which contains the executable "console scan ip"
5. Run the installation command `python run install` (with the `--no-conda` flag if you already have an ASCENT environment)
6. Create a subdirectory called `inputs/sipsource` and place the input `.sip` segmentation files there. See a dummy example in `examples/3D/mock_3D_nerve.sip`.


# Running

1. Create a configuration file by copying `config/templates/3D.json` to `config/<myfile>.json` (replace `<myfile>` with a name of your choice). The configuration file contains all the parameters for the pipeline. The parameters are described in the README file in the `config` directory.
2. Place the necessary ASCENT configuration files in the following locations:
  - `sim.json`: `config` directory
  - `model.json`: `config/models` directory
  - `sample.json`: `config/samples` directory
  - `run.json`: `config` directory
4. Navigate to the `scripts` directory.
5. Call `python 3D_pipeline_ds.py <myfile>.json` (replace with the configuration file you created)
6. The pipeline will run in a folder with the same name as your config file in `datanew`. To generate thresholds, take the directories from the `datanew/<myfile>/7_ascent` folder after running and place them in the relevant ASCENT directories (under the "3d" branch). Add `"post_java_only":true` to your ASCENT run config.

# Deformation

1. Use a `sample.json` where deform is set to `true` (Note: only supports LivaNova cuff and deform_ratio of 1)
2. The pipeline will exit prior to actually running the deformation simulation. You need to do the following:
    - Look through `predeform.sip` to see that no fascicle jump around. If they do, you must manually correct the segmentation.
    - Manually run the deform COMSOL simulation (`datanew/<myfile>/2_slides/doubledeform`). You may need to adjust the contact penalty parameters to get a successful simulation. Even then, you may not succeed. However, if you simulation gets most of the way there, it may still be usable since the ramping is exponential decreasing.
    - Manually export the final simulation step as postdeform.stl
    - Manually Import into Simpleware predeform file and scroll through to check that everything looks ok.
    - Manually use the code in `scripts/sip_postdeform.py` to generate contour coordinates for the deformed fascicles.
    - Run the pipeline again, it should skip deformation now that the proper files are in place.
    - If there are any topological mismatches a window will pop up and ask you to indicate fascicles which watershed segmentation thinks are separate but are not.
    - Once this completes, the pipeline will run as normal.
Note: Sometimes the mesh from COMSOL will be non-manifold. In these cases, you may need to "fix" the import geometry before exporting contours. However, this will result in a LOT of windows popping up during the watershed. To avoid this, write a short code snippet in the API window to loop through and export all possible contours where the mesh is manifold, skipping the others, then apply the fix operation, and export the remaining contours.

Other note: sometimes you may want to run a new cuff on a deformed nerve without redoing the whole deformation process. In this case, follow all the steps to make a new config as normal, then:
1. Copy the "slides" folder to the relevant new directory (in datanew, e.g., from 2Ldeform to 2Ldeformnewcuff)
2. Make sure 3D pipeline has the following settings: slidegen = True, skipsave = True, rundeform = False
3. Run as normal, the pipeline should use the deformed slides and put the new cuff on the geometry
