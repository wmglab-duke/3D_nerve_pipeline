# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 16:44:25 2021

@author: dpm42
"""
import json

from scanip_api3 import *

sipconfig = App.GetInstance().GetInputValue()

# configdict = r'D:\3D_VNS\config\system/configdict.json'

with open(sipconfig) as f:
    config = json.load(f)

masks = ['i', 'n']

# %% Name masks
doc = App.GetInstance().OpenDocument(config['infile'])

doc.GetGenericMaskByName(config['n_mask']).SetName("n")

doc.GetGenericMaskByName(config['i_mask']).SetName("i")

# %% Process masks
doc.ShrinkWrapData(Doc.TargetAllMasks, 10, 10, 10, 10, 10, 10)

doc.GetGenericMaskByName("i").Activate()

doc.MoveMaskTo(doc.GetActiveGenericMask(), 2)

doc.GetGenericMaskByName("n").Activate()

doc.MoveMaskTo(doc.GetActiveGenericMask(), 1)


for m in masks:
    doc.GetGenericMaskByName(m).Activate()

    doc.ApplyCavityFillFilter()

    doc.ApplyIslandRemovalFilter(500)

# doc.ApplyFillGaps(Doc.MostContactSurface, [doc.GetMaskByName("n"), doc.GetMaskByName("i")], True, 0)

#########

# 2022-07-24 10:53:17 - Mask activation
App.GetDocument().GetGenericMaskByName("i").Activate()

# 2022-07-24 10:53:26 - Mask duplication
App.GetDocument().GetActiveGenericMask().Duplicate()

# 2022-07-24 10:53:26 - Mask activation
App.GetDocument().GetGenericMaskByName("Copy of i").Activate()

# 2022-07-24 10:53:30 - Morphological filter
App.GetDocument().ApplyDilateFilter(Doc.TargetMask, 4, 4, 2, 0.0)

# 2022-07-24 10:53:34 - Mask activation
App.GetDocument().GetGenericMaskByName("n").Activate()

# 2022-07-24 10:53:42 - Voxel Boolean
App.GetDocument().ReplaceMaskUsingBooleanExpression(
    "(\"Copy of i\" OR n)",
    App.GetDocument().GetMaskByName("n"),
    App.GetDocument().GetSliceIndices(Doc.OrientationYZ),
    Doc.OrientationYZ,
)

# 2022-07-24 10:53:44 - Mask activation
App.GetDocument().GetGenericMaskByName("Copy of i").Activate()

# 2022-07-24 10:53:47 - Mask removal
App.GetDocument().RemoveMask(App.GetDocument().GetActiveGenericMask())
##########
doc.GetGenericMaskByName("i").Activate()

# 2022-06-13 12:19:20 - Smart mask smoothing
doc.ApplySmartMaskSmoothing([doc.GetMaskByName("i"), doc.GetGenericMaskByName("n")], Doc.Greyscale, 200, True)

doc.ResampleDataByPixelSpacing(
    config['um_per_px'] / 1000 / config['resample_xy'],
    config['um_per_px'] / 1000 / config['resample_xy'],
    config['um_per_slice'] / 1000 / config['resample_factor'],
    Doc.LinearInterpolation,
    Doc.LinearInterpolation,
)


# 2022-06-13 12:32:23 - Morphological filter
App.GetDocument().ApplyOpenFilter(Doc.TargetMask, 2, 2, 0, 0.0)

doc.ApplyFillGaps(Doc.MostContactSurface, [doc.GetMaskByName("n"), doc.GetMaskByName("i")], True, 0)

doc.SaveAs(config['outpath'] + '/debug_prebinary.sip')

# %% Apply binarisation, copy to background, and export
for m in masks:
    doc.GetGenericMaskByName(m).Activate()

    doc.ApplyBinarisationFilter()

    doc.CopyMaskToBackground()

    doc.GetBackgroundByName(m).ImageStackExport(
        config['outpath'] + f"/{m}/", "", "TIFF", Doc.Uint8, Background.Z, False
    )

# 2021-10-26 14:00:20 - Project save
doc.SaveAs(config['outpath'] + '/debug_postbinary.sip')
