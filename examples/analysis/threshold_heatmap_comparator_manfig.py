#!/usr/bin/env python3.7

"""Generate a heatmap of activation thresholds.

The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent.

Note: if more than one heatmap is desired, you must use a Seaborn FacetGrid.
RUN THIS FROM REPOSITORY ROOT
"""
import os

os.chdir("../..")
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from src.core.plotter import heatmaps
from src.core.query import Query

mpl.rcParams["figure.dpi"] = 400
samp2ds = [372, 3729, 3721]
model = 0
simint = 333
samp3ds = [373, 373, 3731]
pairnames = ["Undeformed", "ASCENT", "Structural"]
slidenames = ["3RDS5", "3RdefDS5", "3RdefDS5"]

import matplotlib.colors as mplcolors
import pandas as pd
import seaborn as sns
from matplotlib import colormaps

from src.utils import Object


def add_act_colorbar(ax):
    cmap = colormaps["plasma"]
    cmap.set_bad(color="w")
    cmap = cmap.reversed()
    mappable = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=mplcolors.Normalize(vmin=0, vmax=1),
    )
    cb = plt.colorbar(mappable=mappable, ax=ax, ticks=[0, 1])
    cb.ax.set_ylabel("Percent Activated")


def add_thresh_colorbar(ax, mint, maxt):
    cmap = colormaps["viridis"]
    cmap.set_bad(color="w")
    cmap = cmap.reversed()
    mappable = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=mplcolors.Normalize(vmin=mint, vmax=maxt),
    )
    cb = plt.colorbar(mappable=mappable, ax=ax, ticks=[mint, maxt])
    cb.ax.set_ylabel("Threshold (mA)")


def addpwfd(data, sim):
    with open("examples/analysis/plotconfig.json") as f:
        config = json.load(f)
    nsim_key = config["sim_data"][sim]["nsim_key"]
    for nsim in nsim_key:
        pulse_width = nsim_key[nsim]["pulse_width"]
        fiber_diam = nsim_key[nsim]["fiber_diam"]
        nsim = int(nsim)
        data.loc[data["nsim"] == nsim, "pulse_width"] = pulse_width
        data.loc[data["nsim"] == nsim, "fiber_diam"] = fiber_diam
    return data


plasmap = colormaps["plasma"]
plasmap.set_bad(color="w")
plasmap = plasmap.reversed()
# %%
for samp2d, samp3d, pairname, slidename in zip(samp2ds, samp3ds, pairnames, slidenames):
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
    threshdat = pd.concat([dat2d, dat3d])
    threshdat.query("nsim in [0,5]", inplace=True)
    threshdat = addpwfd(threshdat, "3")
    threshdat.query("nsim==0", inplace=True)
    # %%
    # calculate activation order
    drdat = threshdat.copy().reset_index()  # change this to repeated?
    drdat["percent_activated"] = 0
    drdat = drdat.rename(columns={"sample": "samplenum"})
    for i, row in drdat.iterrows():
        thisdat = drdat.query(
            "threed == @row.threed and samplenum == @row.samplenum and fiber_diam == @row.fiber_diam and sim == @row.sim"
        )
        # percent is number of thresholds less than or equal to this threshold divided by total number of thresholds
        drdat.loc[i, "percent_activated"] = len(thisdat.query("threshold <= @row.threshold")) / len(thisdat)
    drdat = drdat.rename(columns={"samplenum": "sample"})

    # %%plot heatmaps
    titles = [f"2DEM\n{pairname}", "3DM"]
    g = sns.FacetGrid(
        threshdat,
        row="fiber_diam",
        col="sample",
        sharex=False,
        sharey=False,
        margin_titles=True,
    )
    g.map(
        heatmaps,
        *threshdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={"s": 10},
        min_max_ticks=True,
        min_thresh=threshdat.threshold.min(),
        max_thresh=threshdat.threshold.max(),
        colorbar=False,
    )
    minsave = threshdat.threshold.min()
    maxsave = threshdat.threshold.max()
    g.set_xlabels("")
    g.set_ylabels("")
    g.set_titles(col_template="{col_name}", row_template="")
    for ax in g.axes[0, :]:
        thename = titles[0] if ax.get_title()[2] != "3" else titles[1]
        ax.set_title(thename)
    add_thresh_colorbar(
        g.axes,
        threshdat.threshold.min(),
        threshdat.threshold.max(),
    )
    # load contact coordinates
    dire = rf"D:\threed_final\input\contact_coords\{slidename}"
    for i in [1]:
        contact_coords = np.loadtxt(dire + r"\pcs" + str(i + 1) + ".txt", skiprows=8)[:, :2]
        for ax in g.axes.ravel():
            ax.scatter(
                contact_coords[:, 0],
                contact_coords[:, 1],
                s=1,
                color="k",
                label="contacts" if i == 0 else "_",
            )

    # plot activation order
    drdat["threshold"] = drdat["percent_activated"]
    g = sns.FacetGrid(
        drdat,
        row="fiber_diam",
        col="sample",
        sharex=False,
        sharey=False,
        margin_titles=True,
    )
    g.map(
        heatmaps,
        *drdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={"s": 10},
        cmap=plasmap,
        colorbar=False,
    )
    add_act_colorbar(g.axes)
    g.set_xlabels("")
    g.set_ylabels("")
    g.set_titles(col_template="{col_name}", row_template="")
    for ax in g.axes[0, :]:
        thename = titles[0] if ax.get_title()[2] != "3" else titles[1]
        ax.set_title(thename)
    # load contact coordinates
    dire = rf"D:\threed_final\input\contact_coords\{slidename}"
    for i in [1]:
        contact_coords = np.loadtxt(dire + r"\pcs" + str(i + 1) + ".txt", skiprows=8)[:, :2]
        for ax in g.axes.ravel():
            ax.scatter(
                contact_coords[:, 0],
                contact_coords[:, 1],
                s=1,
                color="k",
                label="contacts" if i == 0 else "_",
            )
# %%

samp2d = samp2ds[1]
model = 0
simint = 3
samp3d = samp3ds[1]
titles = [f"2DEM\n{pairnames[1]}", "3DM Undeformed"]

import pandas as pd
import seaborn as sns

from src.utils import Object

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
threshdat = pd.concat([dat2d, dat3d])
# %%
threshdat.query("nsim==0", inplace=True)
g = sns.FacetGrid(threshdat, col="sample", sharex=False, sharey=False)
g.map(
    heatmaps,
    *threshdat.columns,
    sample_object=sample_obj,
    sim_object=sim_obj,
    scatter_kws={"s": 10},
    min_thresh=minsave,
    max_thresh=maxsave,
    colorbar=False,
)
g.set_ylabels("")
g.set_xlabels("")

g.set_titles(col_template="{col_name}", row_template="")
for ax in g.axes[0, :]:
    thename = titles[0] if ax.get_title()[2] != "3" else titles[1]
    ax.set_title(thename)
add_thresh_colorbar(
    g.axes,
    minsave,
    maxsave,
)
# load contact coordinates
dire = rf"D:\threed_final\input\contact_coords\{slidename}"
for i in [1]:
    contact_coords = np.loadtxt(dire + r"\pcs" + str(i + 1) + ".txt", skiprows=8)[:, :2]
    for ax in g.axes.ravel():
        ax.scatter(
            contact_coords[:, 0],
            contact_coords[:, 1],
            s=1,
            color="k",
            label="contacts" if i == 0 else "_",
        )
