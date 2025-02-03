"""Created on Tue Sep  7 16:44:25 2021.

@author: dpm42
"""

import json
import math
import os

from scanip_api3 import *

# temp


# end temp
sipconfig = App.GetInstance().GetInputValue()

sipconfig = r'D:\work\threed\config\system/sipmeshconfig.json' if sipconfig == "" else sipconfig

print(sipconfig)

with open(sipconfig) as f:
    config = json.load(f)

masks = ['i', 'p', 'n']

matpriority = {
    1: 'endoneurium',
    2: 'perineurium',
    3: 'epineurium',
    4: "conductor",
    5: "recess",
    6: "insulator",
    7: "fill",
    8: "medium",
}

# %% import and generate masks
doc = App.GetInstance().ImportStackOfImages(
    App.GetInstance().SearchForImages(config['n_imgs']),
    config['um_per_px'] / 1000,
    config['um_per_px'] / 1000,
    config['um_per_slice'] / 1000,
    CommonImportConstraints().SetWindowLevel(0.0, 0.0),
)

doc.ImportBackgroundFromStackOfImages(
    App.GetInstance().SearchForImages(config['p_imgs']),
    config['um_per_px'] / 1000,
    config['um_per_px'] / 1000,
    config['um_per_slice'] / 1000,
    CommonImportConstraints().SetWindowLevel(0.0, 0.0),
)

doc.ImportBackgroundFromStackOfImages(
    App.GetInstance().SearchForImages(config['i_imgs']),
    config['um_per_px'] / 1000,
    config['um_per_px'] / 1000,
    config['um_per_slice'] / 1000,
    CommonImportConstraints().SetWindowLevel(0.0, 0.0),
)

doc.ImportBackgroundFromStackOfImages(
    App.GetInstance().SearchForImages(config['pcap_imgs']),
    config['um_per_px'] / 1000,
    config['um_per_px'] / 1000,
    config['um_per_slice'] / 1000,
    CommonImportConstraints().SetWindowLevel(0.0, 0.0),
)

doc.GetBackgroundByName("Stack").SetName('n')

doc.GetBackgroundByName("Stack (2)").SetName('p')

doc.GetBackgroundByName("Stack (3)").SetName('i')

doc.GetBackgroundByName("Stack (4)").SetName('pcap')

# 2022-02-16 16:47:38 - Combine backgrounds
doc.ReplaceBackgroundUsingCombineBackgroundsOperation(
    doc.GetBackgroundByName("p"), [doc.GetBackgroundByName("pcap")], doc.GetBackgroundByName("p"), Doc.Maximum
)

for mask in masks:
    doc.GetBackgroundByName(mask).Activate()

    App.GetDocument().CopyBackgroundToMask()

# 2021-10-15 15:24:22 - Flip
# doc.FlipData(Doc.AxisY) nope

for position, mask in enumerate(masks):
    # 2021-10-15 15:25:38 - Mask activation
    doc.GetGenericMaskByName(mask).Activate()

    # 2021-10-15 15:25:38 - Movement of mask
    doc.MoveMaskTo(doc.GetActiveGenericMask(), len(masks) - position)

# 2021-10-15 15:27:55 - Mask activation
doc.GetGenericMaskByName("n").Activate()

# 2021-10-15 15:27:56 - Island removal filter
doc.ApplyIslandRemovalFilter(1000)

# 2021-10-15 15:29:38 - Fill gaps
doc.ApplyFillGaps(Doc.MostContactSurface, [doc.GetMaskByName("p"), doc.GetMaskByName("i")], True, 1000)

# fix for peri too close to epi after end caps
sep_nerve = config['sep_nerve']

pixelval = math.floor(config['sep_nerve'] / config['um_per_px'])

partialval = config['sep_nerve'] / config['um_per_px'] - pixelval

# 2022-05-18 11:27:33 - Mask activation
App.GetDocument().GetGenericMaskByName("p").Activate()

# 2022-05-18 11:27:45 - Mask duplication
App.GetDocument().GetActiveGenericMask().Duplicate()

# 2022-05-18 11:27:46 - Mask activation
App.GetDocument().GetGenericMaskByName("Copy of p").Activate()

# 2022-05-18 11:28:22 - Morphological filter
App.GetDocument().ApplyDilateFilter(Doc.TargetMask, pixelval, pixelval, pixelval, partialval)

# 2022-05-18 11:31:02 - Voxel Boolean
App.GetDocument().ReplaceMaskUsingBooleanExpression(
    "(n OR \"Copy of p\")",
    App.GetDocument().GetMaskByName("n"),
    App.GetDocument().GetSliceIndices(Doc.OrientationXY),
    Doc.OrientationXY,
)

# 2022-05-18 11:32:04 - Mask activation
App.GetDocument().GetGenericMaskByName("Copy of p").Activate()

# 2022-05-18 11:32:05 - Mask removal
App.GetDocument().RemoveMask(App.GetDocument().GetMaskByName("Copy of p"))

# Get the dimensions of the image volume
dims = doc.GetDimensions()
# Store the image-to-global transformation matrix
i2g_mat = doc.GetImageToGlobalTransformationMatrix()
# Define a temporary matrix defining the image-to-new-origin
temp_mat = Matrix.FromTranslation(-dims.GetPhysicalSizeX() / 2, -dims.GetPhysicalSizeY() / 2, 0)
# Temporarily set the image-to-global, generate and export the mesh
doc.SetImageToGlobalTransformationMatrix(temp_mat)
# 2021-09-10 12:57:40 - Import STL file

path = config["stl_path"]
with open(path + '/matmap.json') as f:
    matmap = json.load(f)

matmap["i"] = "endoneurium"
matmap["p"] = "perineurium"
matmap["n"] = "epineurium"

# set material priority
mat_ind = len(matmap)
for key in matpriority:
    for k, v in matmap.items():
        if v == matpriority[key]:
            matmap[k] = {"position": mat_ind, "material": v}
            mat_ind -= 1

stl_files = [x for x in os.listdir(path) if x.endswith('.stl')]
fromstl = [key for key in matmap.keys() if matmap[key]['material'] in ['fill', 'medium']]

# get surfaces
use_nastran = config["mesh"].get("use_nastran", False)
if use_nastran:
    # 2023-05-11 12:12:16 - Import volume mesh file
    App.GetDocument().ImportVolumeMeshFromFile(path + "/alldomain.nas", 0.001, False)

    # 2023-05-11 12:12:24 - Generation of surface(s) from element set(s)
    App.GetDocument().CopyElementSetsToSurfaces(
        App.GetDocument().GetVolumeMeshByName("alldomain").GetElementSets(), True
    )

    # 2023-05-11 12:16:00 - Element set visibility modification
    App.GetDocument().ToggleElementSetVisibility(App.GetDocument().GetVolumeMeshByName("alldomain").GetElementSets())

    # 2023-05-11 12:17:00 - Volume mesh removal
    App.GetDocument().RemoveVolumeMesh(App.GetDocument().GetVolumeMeshByName("alldomain"))

    # get dommap
    with open(path + '/dommap.json') as f:
        dommap = json.load(f)

    # Rename elements to correct surface names
    for setname, setnums in dommap.items():
        for setnum in setnums:
            App.GetDocument().GetSurfaceByName(f"Element set {setnum} (from mesh)").SetName(setname)

    # 2023-05-15 10:48:28 - Surface object removal (removes medium)
    App.GetDocument().RemoveSurface(App.GetDocument().GetSurfaceByName("Element set 1 (from mesh)"))

    for name in fromstl:
        try:
            # 2023-05-15 10:48:28 - Surface object removal (removes fill)
            App.GetDocument().RemoveSurface(App.GetDocument().GetSurfaceByName(name))
        except:
            pass
        doc.ImportSurfaceFromStlFile(path + '/' + name + '.stl', True, 0.001, False)


else:
    for file in stl_files:
        doc.ImportSurfaceFromStlFile(path + '/' + file, True, 0.001, False)


# %% new fill code
# 2021-12-11 13:30:06 - Pad
doc.PadData(1000, 1000, 1000, 1000, 0, 0)

m2s = [key for key in matmap.keys() if matmap[key]['material'] == 'fill']

# 2021-12-11 13:30:20 - Generation of mask(s) from surface object(s)
doc.CopySurfacesToMasks([doc.GetSurfaceByName(name) for name in m2s], doc.AccurateManifold, False)

doc.ShrinkWrapData(Doc.TargetAllMasks, 10, 10, 10, 10, 10, 10)

# 2021-12-11 13:30:20 - Generation of mask(s) from surface object(s)
doc.RemoveSurfaces([doc.GetSurfaceByName(name) for name in m2s], False)

for name in m2s:
    doc.GetGenericMaskByName(name + " (from surface)").SetName(name)
# %% end fill code
# %% begin remesh code

# 2022-05-26 10:01:04 - Remesh surface object
doc.GetSurfaceByName("medium").Remesh(0.5)

surfaces = [os.path.splitext(x)[0] for x in stl_files if os.path.splitext(x)[0] not in m2s]

remsurfs = [s for s in surfaces if s != "medium"]

print(remsurfs)

# 2022-06-09 17:10:51 - Surface model creation
remesh_model = App.GetDocument().CreateSurfaceModel("SurfModel")

# 2022-06-09 17:10:58 - Objects mode activation
App.GetDocument().EnableObjectsMode()

# 2022-06-09 17:11:02 - Part creation
remesh_model.AddSurfaces([App.GetDocument().GetSurfaceByName(x) for x in remsurfs])

# 2022-06-09 17:11:03 - Models mode activation
App.GetDocument().EnableModelsMode()

# 2022-06-09 17:11:27 - Model activation
App.GetDocument().SetActiveModel(remesh_model)

# 2022-06-09 17:11:42 - Model configuration modification
remesh_model.SetExportType(Model.StlFeCfdCad)

# 2022-06-09 17:12:12 - Surface generation
App.GetDocument().GenerateMesh()

# 2022-06-09 17:15:06 - Conversion of meshed part(s) to surface object(s)
remesh_model.CreateSurfacesFromParts([remesh_model.GetPartByName(x) for x in remsurfs])

# 2022-06-09 17:15:09 - Objects mode activation
App.GetDocument().EnableObjectsMode()
for surfname in remsurfs:
    # 2022-06-09 17:15:21 - Rename surface object
    App.GetDocument().GetSurfaceByName(surfname).SetName(surfname + "_imported")

    # 2022-06-09 17:15:27 - Rename surface object
    App.GetDocument().GetSurfaceByName(surfname + " (from mesh)").SetName(surfname)

# %% end remesh code
# 2021-10-15 15:35:29 - FE model creation
femod = doc.CreateFeModel("Model 1")

# 2021-10-15 15:35:35 - Objects mode activation
doc.EnableObjectsMode()

# 2021-10-15 15:35:38 - Part creation
doc.GetModelByName("Model 1").AddMasks(
    [doc.GetGenericMaskByName("n"), doc.GetGenericMaskByName("p"), doc.GetGenericMaskByName("i")]
)


# 2021-12-11 13:48:26 - Part creation
doc.GetModelByName("Model 1").AddMasks([doc.GetGenericMaskByName(name) for name in m2s])

# 2021-10-15 15:35:38 - Models mode activation
doc.EnableModelsMode()

# 2021-10-15 15:35:40 - Objects mode activation
doc.EnableObjectsMode()

# 2021-10-15 15:35:49 - Part creation
femod.AddSurfaces([doc.GetSurfaceByName(surf) for surf in surfaces])

# 2021-10-15 15:35:49 - Models mode activation
doc.EnableModelsMode()

# 2021-10-15 15:36:27 - Model configuration modification
femod.SetUseSmartMaskSmoothing(True)
# 2022-07-15 17:47:53 - Model configuration modification

femod.SetSnapSurfacePartsToModelBounds(True)

# 2022-07-15 17:47:58 - Model configuration modification\
femod.SetNumSmartMaskSmoothingIterations(100)

# 2021-10-15 15:36:30 - Model configuration modification
femod.SetExportUnits(Model.MicronsUnits)

# 2021-10-15 15:37:45 - Contact to boundary creation
femod.AddSurfaceContact(femod.GetPartByName("i"), femod.GetPartByName("medium"))

# 2021-10-15 15:38:14 - Model configuration modification
femod.SetExportType(Model.ComsolNasVolume)

# temp block
# 2022-07-15 17:57:09 - Model configuration modification
femod.SetUseSmallestElementImprovement(True)
# 2022-07-15 18:50:44 - Model configuration modification
femod.SetAdditionalMeshQualityImprovementMaximumOffSurfaceDistance(0.01)
# 2022-07-15 17:57:29 - Model configuration modification
femod.SetSmallestElementImprovementCharacteristicLengthTarget(config['mesh']["max_edge_length"])
# 2022-07-15 17:57:31 - Model configuration modification
femod.SetUseLimitMaximumDisplacement(True)
# 2022-07-15 17:57:37 - Model configuration modification
femod.SetMaximumDisplacementRatio(0.2)
# end temp block

for ob, ob_info in matmap.items():
    # 2021-10-15 15:38:45 - Model configuration modification
    femod.SetEditAdvancedParametersManuallyOnPart(femod.GetPartByName(ob), True)

    if ob_info["material"] == "conductor":
        # 2021-10-15 15:38:52 - Model configuration modification
        femod.SetTargetMinimumEdgeLengthOnPart(femod.GetPartByName(ob), 0.0005)
    else:
        # 2021-10-15 15:38:52 - Model configuration modification
        femod.SetTargetMinimumEdgeLengthOnPart(femod.GetPartByName(ob), config['mesh']['min_edge_length'])

    # 2021-10-15 15:38:55 - Model configuration modification`
    femod.SetTargetMaximumErrorOnPart(femod.GetPartByName(ob), config['mesh']['max_error'])

    # 2021-10-15 15:39:00 - Model configuration modification
    if ob_info["material"] == "conductor":
        femod.SetMaximumEdgeLengthOnPart(femod.GetPartByName(ob), 0.01)
    else:
        femod.SetMaximumEdgeLengthOnPart(femod.GetPartByName(ob), config['mesh']["max_edge_length"])

    # 2021-10-15 15:39:13 - Model configuration modification
    femod.SetInternalChangeRateOnPart(femod.GetPartByName(ob), config['mesh']["internal_change_rate"])

    femod.SetSurfaceChangeRateOnPart(doc.GetActiveModel().GetPartByName(ob), config['mesh']["surface_change_rate"])

    femod.SetTargetNumberElementsAcrossLayerOnPart(
        doc.GetActiveModel().GetPartByName(ob), config['mesh']["n_layer_elements"]
    )

    femod.SetSmoothAgainstBackground(config['mesh']["smooth_against_background"])

    try:
        femod.GetPartsContainer().GetPartByName(ob).MoveTo(ob_info["position"])
    except Exception:
        pass

    femod.GetPartByName(ob).SetMaterial(PlaceholderMaterial(ob_info["material"]))

# 2022-04-26 16:26:12 - Model configuration modification
if config['mesh'].get('second_order') is True:
    femod.SetHigherOrder(True)
    femod.SetCurvedEdges(True)
else:
    print('Warning: Using first order elements')

# required twice
for ob, ob_info in matmap.items():
    try:
        femod.GetPartsContainer().GetPartByName(ob).MoveTo(ob_info["position"])
    except Exception:
        pass

doc.GetGenericMaskByName("p").Activate()

# 2022-06-13 12:32:23 - Morphological filter
App.GetDocument().ApplyOpenFilter(Doc.TargetMask, 2, 2, 0, 0.0)

# 2021-11-02 12:40:33 - Isolated cavity and island removal
doc.RemoveIsolatedCavitiesAndIslands([doc.GetMaskByName("n"), doc.GetMaskByName("p"), doc.GetMaskByName("i")], 101, 101)

doc.ShrinkWrapData(Doc.TargetAllMasks, 10, 10, 10, 10, 10, 10)

doc.SaveAs(config['outpath'] + '/mesh_debug.sip')
if config["run_type"] == "cluster":
    # 2021-10-26 14:00:20 - Project save
    doc.GenerateMesh()
    # 2021-12-17 08:44:39 - COMSOL export
    doc.ExportComsolNasVolume(config["outpath"] + '/mesh.nas', False)
    # 2021-10-26 14:00:20 - Project save
    doc.SaveAs(config['outpath'] + '/mesh.sip')
