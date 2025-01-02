#!/usr/bin/env python3.7

"""Generate a heatmap of activation thresholds.

The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent.

Note: if more than one heatmap is desired, you must use a Seaborn FacetGrid.
RUN THIS FROM REPOSITORY ROOT
"""
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

os.chdir("../..")

from src.core.plotter import heatmaps
from src.core.query import Query
from src.utils import Object

samp2d = 5701
model = 0
simint = 330
samp3d = 5731


q = Query(
    {
        "partial_matches": True,
        "include_downstream": True,
        "indices": {"sample": [samp2d], "model": [model], "sim": [simint]},
    }
).run()
dat2d = q.data(thresh_only=True)
dat2d["threed"] = False
q3 = Query(
    {
        "partial_matches": True,
        "include_downstream": True,
        "indices": {"sample": [samp3d], "model": [model], "sim": [simint]},
    }
).run()
dat3d = q3.data(source_sample=samp2d, thresh_only=True)
dat3d["threed"] = True
sample_obj = q.get_object(Object.SAMPLE, [samp2d])
sim_obj = q.get_object(Object.SIMULATION, [samp2d, model, simint])
threshdat = pd.concat([dat2d, dat3d]).query("fiberset_index==1")
# %%
g = sns.FacetGrid(threshdat, row="active_src_index", col="sample", sharex=False, sharey=False, margin_titles=True)
g.map(heatmaps, *threshdat.columns, sample_object=sample_obj, sim_object=sim_obj, scatter_kws={"s": 25})

save_directory = os.path.join("output", "analysis")
os.makedirs(save_directory, exist_ok=True)
plt.savefig(os.path.join(save_directory, "threshold_heatmap.png"), dpi=400, bbox_inches="tight")
