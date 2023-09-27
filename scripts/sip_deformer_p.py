# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 16:44:25 2021

@author: dpm42
Notes:
after this test goes test whether exhaustive gives better results with nerve-nerve
then if works better add to code
also need to perhaps discard no shifting in z assumption
need to find a way to add inelastic (check if modulus is still changed)
Can also try again with nerve-nerve separation
Note also compare all parameters to the manual example I did
Apply more squeeze to nerve than sep so then there is wiggle room for separation of fascicles
can also grow fascicles before doing the sep and then shrink them after for separation
May need to change contact type to symmetric
"""
import json
import os

from scanip_api3 import *

sipconfig = App.GetInstance().GetInputValue()

sipconfig = r'D:\work\threed\config\system\sipdefconfig.json' if sipconfig == "" else sipconfig

with open(sipconfig) as f:
    config = json.load(f)

# %% Load doc
# 2022-12-15 14:48:23 - Import stack of images
# Please note, if wanting to import all image files in the directory
# (without cropping) this more compact code may be used instead:
App.GetInstance().ImportStackOfImages(
    App.GetInstance().SearchForImages(os.path.join(config['out'], "p")),
    0.005,
    0.005,
    0.02,
    CommonImportConstraints().SetWindowLevel(0.0, 0.0),
)  # TODO: make this not hard coded

doc = App.GetDocument()

# 2022-12-15 14:48:33 - Background to mask copy
doc.CopyBackgroundToMask()

# 2022-12-15 14:48:46 - Mask renaming
doc.GetGenericMaskByName("Stack").SetName("p")

App.GetInstance().GetActiveDocument().ImportBackgroundFromStackOfImages(
    App.GetInstance().SearchForImages(os.path.join(config['out'], "n")),
    0.005,
    0.005,
    0.02,
    CommonImportConstraints().SetWindowLevel(0.0, 0.0),
)
# 2022-12-15 14:48:33 - Background to mask copy
doc.CopyBackgroundToMask()

# 2022-12-15 14:48:46 - Mask renaming
doc.GetGenericMaskByName("Stack (2)").SetName("n")

os.makedirs(config['out'], exist_ok=True)
#%% Setup
zmove = config['transition_distance'] / doc.GetDimensions().GetSpacingZ() / 1000

# crop
doc.CropData(
    0,
    doc.GetDimensions().GetVoxelCountX(),
    0,
    doc.GetDimensions().GetVoxelCountY(),
    int(config['span'][0] - zmove),
    int(config['span'][1] + zmove),
)

# 2022-12-21 12:37:18 - Mask activation
App.GetDocument().GetGenericMaskByName("p").Activate()

# # 2022-12-21 12:37:53 - Morphological filter
# App.GetDocument().ApplyOpenFilter(Doc.TargetMask, 3, 3, 3, 0.0)

# # 2022-12-21 12:38:05 - Morphological filter
# App.GetDocument().ApplyOpenFilter(Doc.TargetMask, 3, 3, 0, 0.0)

#%% mesh generation
smod = doc.CreateSurfaceModel("Model 1")

smod.AddMask(App.GetDocument().GetGenericMaskByName("p"))

smod.SetUseSmartMaskSmoothing(True)

smod.SetNumSmartMaskSmoothingIterations(100)

smod.SetExportType(Model.StlFeCfdCad)

smod.SetCompoundCoarsenessOnPart(smod.GetPartByName("p"), -50)

#%% mesh export
dims = doc.GetDimensions()
# Store the image-to-global transformation matrix
i2g_mat = doc.GetImageToGlobalTransformationMatrix()
# Define a temporary matrix defining the image-to-new-origin
temp_mat = Matrix.FromTranslation(-dims.GetPhysicalSizeX() / 2, -dims.GetPhysicalSizeY() / 2, 0)
# Temporarily set the image-to-global, generate and export the mesh
doc.SetImageToGlobalTransformationMatrix(temp_mat)

doc.GenerateMesh()

# save
doc.SaveAs(os.path.join(config['out'], "deformed.sip"))

# 2021-09-10 12:57:40 - Import STL file
doc.ExportStl(os.path.join(config['out'], "predeform.stl"), Doc.WriteSingleBinaryFile, False)
