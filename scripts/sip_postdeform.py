# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 16:44:25 2021

@author: dpm42
"""

import json
import os
import sys

from scanip_api3 import *

sipconfig = App.GetInstance().GetInputValue()

sipconfig = r'D:/work/threed/datanew/5RDS5def/2_slides/doubledeform/sipdefconfig.json' if sipconfig == "" else sipconfig

with open(sipconfig) as f:
    config = json.load(f)

# %% Load doc
doc = App.GetInstance().OpenDocument(os.path.join(config['out'], "deformed.sip"))

# 2022-12-15 18:21:50 - Import STL file
doc.ImportMultipleSurfacesFromStlFile(os.path.join(config['out'], "postdeform.stl"), False, 1.0, False)

# 2022-12-15 18:22:01 - Surface object removal
App.GetDocument().RemoveSurface(App.GetDocument().GetSurfaceByName("postdeform_part1"))

# 2022-12-15 18:22:01 - Surface object removal
App.GetDocument().RemoveSurface(App.GetDocument().GetSurfaceByName("postdeform_part2"))

# 2023-01-21 17:20:43 - Objects mode activation
App.GetDocument().EnableObjectsMode()

# 2023-01-21 17:21:10 - Mask activation
App.GetDocument().GetGenericMaskByName("n").Activate()

# 2023-01-21 17:21:14 - Surface model creation
App.GetDocument().CreateSurfaceModel("Model 2")

# 2023-01-21 17:21:16 - Objects mode activation
App.GetDocument().EnableObjectsMode()

# 2023-01-21 17:21:16 - Part creation
App.GetDocument().GetModelByName("Model 2").AddMask(App.GetDocument().GetGenericMaskByName("n"))

# 2023-01-21 17:21:17 - Models mode activation
App.GetDocument().EnableModelsMode()

# 2023-01-21 17:21:23 - Fast preview
App.GetDocument().GenerateFastPreview()

# 2023-01-21 17:21:42 - Surface generation
App.GetDocument().GenerateMesh()

# 2023-01-21 17:21:51 - Conversion of meshed part(s) to surface object(s)
App.GetDocument().GetActiveModel().CreateSurfacesFromParts([App.GetDocument().GetActiveModel().GetPartByName("n")])

# 2023-01-21 17:21:53 - Objects mode activation
App.GetDocument().EnableObjectsMode()

# 2023-01-21 17:21:54 - Surface object activation
App.GetDocument().GetSurfaceByName("n (from mesh)").Activate()

# 2023-01-21 17:21:54 - Surface object activation
App.GetDocument().GetSurfaceByName("postdeform_part3").Activate()

doc.SaveAs(config['outpath'] + '/postdeformedprecontour.sip')

sys.exit('Perform the rest manually for now')

import os

from scanip_api3 import *

# TODO: set x-y plane as active
# todo: delete dirs before export
# TODO: currently need to manually open file and run on a separate directory for the nerve and fascicles

# Obtain a reference to the application
app = App.GetInstance()

# Display the ScanIP version number in a message box
app.ShowMessage(app.GetVersion(), 'ScanIP version')

dire = r'D:\work\threed\datanew\5RDS5def\2_slides\doubledeform\postdeformtxts\n'

view = app.GetDocument().GetActiveSliceView()

contours = app.GetDocument().GetCrossSectionContours()

for i in range(app.GetDocument().GetDimensions().GetVoxelCountZ()):
    view.SetActiveSlice(i)
    os.makedirs(os.path.join(dire, f'{i}'), exist_ok=True)

    contours.CreateFromSliceView(view)

    allc = contours.GetAll()
    for c in allc:
        pts = c.GetControlPoints()
        vertices_list = [[pt.GetX(), pt.GetY()] for pt in pts]
        # save vertices_list to file
        with open(os.path.join(dire, f'{i}', f'{c.GetName()}.txt'), 'w') as f:
            f.write(str(vertices_list))
        contours.DeleteByName(c.GetName())


app.ShowMessage('done')

import os

from scanip_api3 import *

# TODO: set x-y plane as active
# TODO: currently need to manually open file and run on a separate directory for the nerve and fascicles

# Obtain a reference to the application
app = App.GetInstance()

# Display the ScanIP version number in a message box
app.ShowMessage(app.GetVersion(), 'ScanIP version')

dire = r'D:\work\threed\datanew\<insert>\2_slides\doubledeform\postdeformtxts\p'

view = app.GetDocument().GetActiveSliceView()

contours = app.GetDocument().GetCrossSectionContours()

for i in range(app.GetDocument().GetDimensions().GetVoxelCountZ()):
    view.SetActiveSlice(i)
    os.makedirs(os.path.join(dire, f'{i}'), exist_ok=True)

    contours.CreateFromSliceView(view)

    allc = contours.GetAll()
    for c in allc:
        pts = c.GetControlPoints()
        vertices_list = [[pt.GetX(), pt.GetY()] for pt in pts]
        # save vertices_list to file
        with open(os.path.join(dire, f'{i}', f'{c.GetName()}.txt'), 'w') as f:
            f.write(str(vertices_list))
        contours.DeleteByName(c.GetName())

# doc.SaveAs(config['outpath'] + '/postdeformed.sip')
