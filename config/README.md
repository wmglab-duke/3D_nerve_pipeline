# Configuration Documentation

This document describes the available config parameters for the 3D pipeline. Parameters with TODO have immediate action items.

## General Configuration

- **sourcedir**: Directory relative to the current working directory where source files are located. Example: `"../sipsource"` indicates a folder one level up in the directory hierarchy. #TODO remove in favor of assuming this value

- **sourcefile**: Name of the simpleware source file for the process. Example: `"2L_strip.sip"`.

- **path**: Directories for different stages or aspects of the workflow, categorized as follows: #TODO remove this and migrate to ASCENT data structure
  - `slides`: Directory for slide images.
  - `fibers`: Directory for fiber datasets.
  - `mesh`: Directory for meshing results.
  - `comsol`: Directory for COMSOL simulation files.
  - `ascent`: Directory for ascent data.

- **config**: Path and names of essential configuration files: #TODO remove this and migrate to ASCENT data structure
  - `path`: Directory containing configuration files.
  - `run`: Run parameters file.
  - `model`: Model specifications file.
  - `sample`: Sample information file.

- **project_path**: Absolute path to the project's base directory on the system. Example: `"D:\\threed_final\\datanew\\2Lstrip"`. #TODO remove in favor of assuming this value

## Input Parameters

- **i_mask**: Endoneurium input mask in Simpleware file (note must be set as visible or an error will occur.)
- **n_mask**: Epineurium input mask in Simpleware file (note must be set as visible or an error will occur.)
- **um_per_px**: Micrometers per pixel for input mask resolution.
- **um_per_slice**: Slice thickness in micrometers for z-resolution in 3D reconstructions.
- **z+ / z-**: Directional indicators for the caudal and rostral ends, respectively. Only for user reference, not used in pipeline.

## Preprocessing Parameters

- **upsample_z**: Factor for upsampling data in the z-direction. For example, `"2"` doubles the z-resolution, while `"0.5"` halves it.
- **buffer**: Additional open space added around the nerve during image export (keeps the nerve from touching the edge of the image). Units of micrometers.
- **crop**: Indicates whether cropping is enabled (`"True"` for active cropping).
- **center_target**: Target structure to center during preprocessing, e.g., `"nerve"`.
- **nsmoothing**: Degree of smoothing for the n_mask.
- **fsmoothing**: Degree of smoothing for the fiber dataset.

## Output Parameters

- **image_um_per_px**: Spatial resolution of output images in micrometers per pixel.
- **id_font**: Font size for identification labels in output images or diagrams.
- **addl_cuff_buffer**: Additional buffer space around cuffs in the output.
- **slices**: Positions in micrometers of key slices or sections (`rostral`, `center`, `caudal`).

## Mesh Generation Parameters

- **min_edge_length**: Minimum edge length in the generated mesh.
- **max_error**: Maximum error tolerance in mesh generation.
- **max_edge_length**: Maximum edge length for controlling mesh density.
- **internal_change_rate**: Rate of change for internal mesh parameters.
- **surface_change_rate**: Rate of change for surface mesh parameters.
- **n_layer_elements**: Number of elements in the normal direction for mesh layering.
- **smooth_against_background**: Boolean indicating if mesh smoothing considers background elements.
- **second_order**: Specifies use of second-order elements for increased accuracy.
- **use_nastran**: Boolean indicating whether to use Nastran for mesh generation or analysis.

This documentation outlines the basic functionality and purpose of each parameter in the configuration. For more detailed explanations or specific use cases, additional context provided by the user would be beneficial.
# Configuration Documentation

This document describes the parameters used in the configuration JSON object for a specific scientific or engineering process involving image processing, mesh generation, and simulations.

## General Configuration #TODO migrate to ASCENT run config

- **sourcedir**: Directory relative to the current working directory where source files are located. Example: `"../sipsource"` indicates a folder one level up in the directory hierarchy.

- **sourcefile**: Name of the primary source file for the process. Example: `"2L_strip.sip"` suggests a specific format or data relevant to the workflow.

- **path**: Directories for different stages or aspects of the workflow, categorized as follows:
  - `slides`: Directory for slide images.
  - `fibers`: Directory for fiber datasets.
  - `mesh`: Directory for meshing results.
  - `comsol`: Directory for COMSOL simulation files.
  - `ascent`: Directory for ascent data.

- **config**: Path and names of essential configuration files:
  - `path`: Directory containing configuration files.
  - `run`: Run parameters file.
  - `model`: Model specifications file.
  - `sample`: Sample information file.

- **project_path**: Absolute path to the project's base directory on the system. Example: `"D:\\threed_final\\datanew\\2Lstrip"`.

## Input Parameters #TODO migrate to ASCENT sample config.

- **i_mask**: Identifier for the endoneurium input mask after opening operations.
- **n_mask**: Identifier for the epineurium input mask after a destepping process.
- **um_per_px**: Micrometers per pixel for input image resolution.
- **um_per_slice**: Slice thickness in micrometers for z-resolution in 3D reconstructions.
- **z+ / z-**: Directional indicators for the caudal and rostral ends, respectively.

## Preprocessing Parameters #TODO migrate to ASCENT sample config

- **upsample_z**: Factor for upsampling data in the z-direction. Note, in the current iteration, the preprocessing script will not run again if output images exist from preprocessing, even if you change this parameter.
- **buffer**: Size of the buffer around the target area during preprocessing.
- **nsmoothing**: Degree of smoothing for the epineurium. #TODO make this come from ASCENT sample config
- **fsmoothing**: Degree of smoothing for the fascicles. #TODO make this come from ASCENT sample config

## Output Parameters #TODO migrate to ASCENT sample config

- **image_um_per_px**: Spatial resolution of output images in micrometers per pixel after preprocessing.
- **id_font**: Font size for identification labels in output images or diagrams.
- **addl_cuff_buffer**: Additional space between cuff and nerve (radius, units of micrometers).
- **slices**: Dictionary defining output slices for use in generating 2D extrusion models. Keys should be the name and values should be the position along the nerve in micrometers.

## Mesh Generation Parameters #TODO migrate to ASCENT model config

- **min_edge_length**: Minimum edge length in the generated mesh. See simpleware documentation for more details.
- **max_error**: Maximum error tolerance in mesh generation. See simpleware documentation for more details.
- **max_edge_length**: Maximum edge length for controlling mesh density. See simpleware documentation for more details.
- **internal_change_rate**: Rate of change for internal mesh parameters. See simpleware documentation for more details.
- **surface_change_rate**: Rate of change for surface mesh parameters. See simpleware documentation for more details.
- **n_layer_elements**: Minimum number of elements across thin layers. See simpleware documentation for more details.
- **smooth_against_background**: Recommend true, makes mask intersections with stls less janky. See simpleware documentation for more details.
- **second_order**: Specifies use of second-order elements for increased accuracy. Recommened true. See simpleware documentation for more details.
- **use_nastran**: Boolean indicating whether to use Nastran for mesh generation or analysis. When exporting cuff geometry, COMSOL will currently export stls and a nastran file. This should be true for most cases, but the LivaNova cuffs do not currently work with this option set to true. Your selection here defines what kind of file simpleware will try to use for importing cuffs.
