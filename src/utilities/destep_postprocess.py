from scanip_api3 import *

# 2022-06-16 14:52:31 - Background to mask copy
App.GetDocument().CopyBackgroundToMask()

# 2022-06-16 14:52:34 - Background removal
App.GetDocument().RemoveBackground(App.GetDocument().GetBackgroundByName("Stack"))

# 2022-06-16 14:52:38 - Mask renaming
App.GetDocument().GetGenericMaskByName("Stack").SetName("endoneurium")

# 2022-06-16 14:52:40 - Mask activation
App.GetDocument().GetGenericMaskByName("endoneurium_destepped").Activate()

# 2022-06-16 14:52:41 - Mask duplication
App.GetDocument().GetActiveGenericMask().Duplicate()

# 2022-06-16 14:52:41 - Mask activation
App.GetDocument().GetGenericMaskByName("Copy of endoneurium_destepped").Activate()

# 2022-06-16 14:52:58 - Mask renaming
App.GetDocument().GetGenericMaskByName("Copy of endoneurium_destepped").SetName("endoneurium_ds_opened")

# 2022-06-16 14:53:22 - Resample of data by pixel spacing (mm)
App.GetDocument().ResampleDataByPixelSpacing(0.01, 0.01, 0.02, Doc.LinearInterpolation, Doc.LinearInterpolation)

# 2022-06-16 14:53:22 - Mask visibility toggling
App.GetDocument().GetActiveGenericMask().SetVisible(True)

# 2022-06-16 14:53:23 - Mask visibility toggling
App.GetDocument().GetGenericMaskByName("endoneurium").SetVisible(False)

# 2022-06-16 14:53:42 - Morphological filter
App.GetDocument().ApplyOpenFilter(Doc.TargetMask, 2, 2, 2, 0.0)

# 2022-06-16 14:53:46 - Mask activation
App.GetDocument().GetGenericMaskByName("endoneurium_destepped").Activate()

# 2022-06-16 14:53:46 - Movement of mask
App.GetDocument().MoveMaskTo(App.GetDocument().GetActiveGenericMask(), 3)

# 2022-06-16 14:53:50 - Mask visibility toggling
App.GetDocument().GetActiveGenericMask().SetVisible(True)

# 2022-06-16 14:53:51 - Mask visibility toggling
App.GetDocument().GetGenericMaskByName("endoneurium").SetVisible(True)

# 2022-06-16 14:53:54 - Mask visibility toggling
App.GetDocument().GetGenericMaskByName("endoneurium").SetVisible(False)

# 2022-06-16 14:53:56 - Mask activation
App.GetDocument().GetGenericMaskByName("endoneurium_ds_opened").Activate()

# 2022-06-16 14:53:57 - Movement of mask
App.GetDocument().MoveMaskTo(App.GetDocument().GetActiveGenericMask(), 3)
