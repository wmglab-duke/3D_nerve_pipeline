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

- **project_path**: Absolute path to the project's base directory on the system. Example: `"D:\\threed_final\\datanew\\2Lstrip"`.

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
