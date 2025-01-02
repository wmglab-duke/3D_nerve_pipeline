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
import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colormaps
from src.core.plotter import datamatch_merge, heatmaps
from src.core.query import Query
from src.utils import Object

# Set matplotlib parameters
mpl.rcParams["figure.dpi"] = 400
sns.set(style='white', context='paper')

# Constants
SAMP2DS = [3721]
SAMP3DS = [3731]
PAIRNAMES = ["Structural"]
SLIDENAMES = ["3RdefDS5"]
MODEL = 0
SIMINT_DEFAULT = 333
CONTACT_DIR_TEMPLATE = r"D:\threed_final\input\contact_coords\{}"
PLOT_CONFIG_PATH = "examples/analysis_daniel/plotconfig.json"


def change_working_directory():
    """Change the current working directory to the repository root."""
    pass


def add_colorbar(ax, cmap_name, label, vmin, vmax, reversed_cmap=True, ticks=None):
    """Add a colorbar to the given axis.

    Parameters:
        ax (matplotlib.axes.Axes): The axis to attach the colorbar to.
        cmap_name (str): Name of the colormap.
        label (str): Label for the colorbar.
        vmin (float): Minimum value for normalization.
        vmax (float): Maximum value for normalization.
        reversed_cmap (bool): Whether to reverse the colormap.
        ticks (list): Tick positions on the colorbar.
    """
    cmap = colormaps[cmap_name].copy()
    cmap.set_bad(color="w")
    if reversed_cmap:
        cmap = cmap.reversed()
    norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)
    mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = plt.colorbar(mappable=mappable, ax=ax, ticks=ticks)
    cb.ax.set_ylabel(label)


def add_act_colorbar(ax):
    """Add a colorbar for activation percentage."""
    add_colorbar(ax=ax, cmap_name="plasma", label="Percent Activated", vmin=0, vmax=1, reversed_cmap=True, ticks=[0, 1])


def add_thresh_colorbar(ax, mint, maxt, cmap='viridis'):
    """Add a colorbar for threshold values."""
    add_colorbar(
        ax=ax, cmap_name="viridis", label="Threshold (mA)", vmin=mint, vmax=maxt, reversed_cmap=True, ticks=[mint, maxt]
    )


def add_pwfd(data, sim):
    """Add pulse width and fiber diameter information to the data.

    Parameters:
        data (pd.DataFrame): The dataframe to update.
        sim (str): Simulation identifier.

    Returns:
        pd.DataFrame: Updated dataframe with 'pulse_width' and 'fiber_diam'.
    """
    with open(PLOT_CONFIG_PATH) as f:
        config = json.load(f)
    nsim_key = config["sim_data"][sim]["nsim_key"]
    for nsim, params in nsim_key.items():
        nsim_int = int(nsim)
        data.loc[data["nsim"] == nsim_int, ["pulse_width", "fiber_diam"]] = params["pulse_width"], params["fiber_diam"]
    return data


def load_contact_coords(slidename, index=2):
    """Load contact coordinates from a specified file.

    Parameters:
        slidename (str): Slide name to construct the file path.
        index (int): Index of the contact file.

    Returns:
        np.ndarray: Array of contact coordinates.
    """
    dire = CONTACT_DIR_TEMPLATE.format(slidename)
    filepath = os.path.join(dire, f"pcs{index}.txt")
    return np.loadtxt(filepath, skiprows=8)[:, :2]


def calculate_percent_activated(threshdat):
    """Calculate the percent of activation for each threshold.

    Parameters:
        threshdat (pd.DataFrame): Data containing threshold information.

    Returns:
        pd.DataFrame: DataFrame with an additional 'percent_activated' column.
    """
    drdat = threshdat.copy().reset_index()
    drdat["percent_activated"] = 0
    drdat = drdat.rename(columns={"sample": "samplenum"})

    for i, row in drdat.iterrows():
        mask = (
            (drdat["threed"] == row.threed)
            & (drdat["samplenum"] == row.samplenum)
            & (drdat["fiber_diam"] == row.fiber_diam)
            & (drdat["sim"] == row.sim)
        )
        relevant = drdat[mask]
        drdat.at[i, "percent_activated"] = (relevant["threshold"] <= row.threshold).sum() / len(relevant)

    return drdat.rename(columns={"samplenum": "sample"})


def plot_heatmaps(threshdat, sample_obj, sim_obj, titles, min_thresh, max_thresh, slidename):
    """Plot heatmaps using Seaborn FacetGrid.

    Parameters:
        threshdat (pd.DataFrame): Data for plotting.
        sample_obj: Sample object from the query.
        sim_obj: Simulation object from the query.
        titles (list): Titles for the heatmaps.
        min_thresh (float): Minimum threshold value.
        max_thresh (float): Maximum threshold value.
        slidename (str): Slide name for loading contact coordinates.
    """
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
        min_thresh=min_thresh,
        max_thresh=max_thresh,
        colorbar=False,
    )

    g.set_xlabels("")
    g.set_ylabels("")
    g.set_titles(col_template="{col_name}", row_template="")

    for ax in g.axes[0, :]:
        # Assuming the title's third character indicates the type
        thename = titles[0] if ax.get_title()[2] != "3" else titles[1]
        ax.set_title(thename)

    add_thresh_colorbar(g.axes, min_thresh, max_thresh)

    # Load and plot contact coordinates
    contact_coords = load_contact_coords(slidename)
    for ax in g.axes.ravel():
        ax.scatter(contact_coords[:, 0], contact_coords[:, 1], s=1, color="k", label="contacts")

    # now on a single plot, plot residual
    # create new dataframe where threshold is 3D threshold - extrusion threshold
    resdat = datamatch_merge(
        threshdat.query('threed==True'),
        threshdat.query('threed==False'),
        "threshold",
        merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
    ).drop(columns="threed")
    resdat["threshold2d"] = resdat["threshold"]
    resdat["threshold"] = resdat["threshold3d"] - resdat["threshold2d"]
    max_abs_thresh = np.percentile(resdat.threshold.abs(), 90)

    g = sns.FacetGrid(
        resdat,
        row="fiber_diam",
        sharex=False,
        sharey=False,
    )
    g.map(
        heatmaps,
        *resdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={"s": 10},
        min_max_ticks=True,
        min_thresh=-max_abs_thresh,
        max_thresh=max_abs_thresh,
        colorbar=True,
        cmap=plt.cm.get_cmap("PiYG"),
    )

    for ax in g.axes.ravel():
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title('Threshold\n(True-3D - Extrusion)')

    # now on a single plot, percent difference
    # create new dataframe where threshold is 3D threshold - extrusion threshold
    resdat["threshold"] = 100 * (resdat["threshold3d"] - resdat["threshold2d"]) / resdat["threshold3d"]
    max_abs_thresh = np.percentile(resdat.threshold.abs(), 90)

    g = sns.FacetGrid(
        resdat,
        row="fiber_diam",
        sharex=False,
        sharey=False,
    )
    g.map(
        heatmaps,
        *resdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={"s": 10},
        min_max_ticks=True,
        min_thresh=-15,
        max_thresh=15,
        colorbar=True,
        cmap=plt.cm.get_cmap("PiYG"),
        cbar_kws={'title': '%'},
    )

    for ax in g.axes.ravel():
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title('Threshold\npd(True-3D, Extrusion)')

    # now go back to the original data, and do a plot where each inner is normalized by its mean
    normdat = threshdat.copy()
    # first need to assign the "inner" value to threed=True by matching it to the threed=False value
    for row in normdat.itertuples():
        if row.threed:
            normdat.at[row.Index, "inner"] = normdat.query(
                f'nsim=={row.nsim} and fiber_diam=={row.fiber_diam} and master_fiber_index=={row.master_fiber_index} and threed==False'
            ).inner.values[0]
    normdat["threshold"] /= normdat.groupby(["sample", "inner", "nsim", "model"])["threshold"].transform("mean")
    g = sns.FacetGrid(
        normdat,
        row="fiber_diam",
        col="sample",
        sharex=False,
        sharey=False,
        margin_titles=True,
    )
    g.map(
        heatmaps,
        *normdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={"s": 10},
        min_max_ticks=True,
        min_thresh=0.8,
        max_thresh=1.2,
        colorbar=True,
        cmap=plt.cm.get_cmap("PiYG"),
        cbar_kws={'title': '', 'label': 'Normalized Threshold'},
    )
    g.set_titles(col_template='', row_template='')
    g.set_xlabels('')
    g.set_ylabels('')

    # add_thresh_colorbar(g.axes, 0.9, 1.1,cmap='PiYG')


def plot_activation_order(drdat, sample_obj, sim_obj, titles, slidename, plasmap):
    """Plot activation order heatmaps.

    Parameters:
        drdat (pd.DataFrame): Data with percent activated.
        sample_obj: Sample object from the query.
        sim_obj: Simulation object from the query.
        titles (list): Titles for the heatmaps.
        slidename (str): Slide name for loading contact coordinates.
        plasmap: Reversed plasma colormap.
    """
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

    # Load and plot contact coordinates
    contact_coords = load_contact_coords(slidename)
    for ax in g.axes.ravel():
        ax.scatter(contact_coords[:, 0], contact_coords[:, 1], s=1, color="k", label="contacts")

    # same as above, plot residual
    resdat = datamatch_merge(
        drdat.query('threed==True'),
        drdat.query('threed==False'),
        "threshold",
        merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
    ).drop(columns="threed")
    resdat["threshold2d"] = resdat["threshold"]
    resdat["threshold"] = resdat["threshold3d"] - resdat["threshold2d"]
    max_abs_thresh = np.percentile(resdat.threshold.abs(), 90)

    g = sns.FacetGrid(
        resdat,
        row="fiber_diam",
        sharex=False,
        sharey=False,
    )
    g.map(
        heatmaps,
        *resdat.columns,
        sample_object=sample_obj,
        sim_object=sim_obj,
        scatter_kws={"s": 10},
        min_max_ticks=True,
        min_thresh=-max_abs_thresh,
        max_thresh=max_abs_thresh,
        colorbar=True,
        cmap=plt.cm.get_cmap("seismic"),
    )

    for ax in g.axes.ravel():
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title('Activation\n(True-3D - Extrusion)')


def process_sample(samp2d, samp3d, pairname, slidename, model, simint):
    """Process a single sample by querying data, preparing it, and plotting heatmaps.

    Parameters:
        samp2d (int): 2D sample identifier.
        samp3d (int): 3D sample identifier.
        pairname (str): Name of the sample pair.
        slidename (str): Slide name for contact coordinates.
        model (int): Model identifier.
        simint (int): Simulation interval identifier.
    """
    # Query 2D data
    q = Query(
        {
            "partial_matches": True,
            "include_downstream": True,
            "indices": {"sample": [samp2d], "model": [model], "sim": [simint]},
        }
    ).run()
    dat2d = q.data(thresh_only=True)
    dat2d["threed"] = False

    # Query 3D data
    q3 = Query(
        {
            "partial_matches": True,
            "include_downstream": True,
            "indices": {"sample": [samp3d], "model": [model], "sim": [simint]},
        }
    ).run()
    dat3d = q3.data(source_sample=samp2d, thresh_only=True)
    dat3d["threed"] = True

    # Combine data
    threshdat = pd.concat([dat2d, dat3d])
    threshdat = threshdat[threshdat["nsim"].isin([0, 5])]
    threshdat = add_pwfd(threshdat, "3")
    threshdat = threshdat[threshdat["nsim"] == 0]

    # Calculate activation percentage
    drdat = calculate_percent_activated(threshdat)

    # Retrieve objects for plotting
    sample_obj = q.get_object(Object.SAMPLE, [samp2d])
    sim_obj = q.get_object(Object.SIMULATION, [samp2d, model, simint])

    # Plot heatmaps
    titles = [f"Extrusion", "True-3D"]
    plot_heatmaps(
        threshdat=threshdat,
        sample_obj=sample_obj,
        sim_obj=sim_obj,
        titles=titles,
        min_thresh=threshdat.threshold.min(),
        max_thresh=threshdat.threshold.max(),
        slidename=slidename,
    )

    # Plot activation order
    plasmap = colormaps["plasma"].reversed()
    plot_activation_order(
        drdat=drdat, sample_obj=sample_obj, sim_obj=sim_obj, titles=titles, slidename=slidename, plasmap=plasmap
    )


def main():
    """Main function to execute the heatmap generation."""
    change_working_directory()

    # Process each sample pair
    for samp2d, samp3d, pairname, slidename in zip(SAMP2DS, SAMP3DS, PAIRNAMES, SLIDENAMES):
        process_sample(
            samp2d=samp2d, samp3d=samp3d, pairname=pairname, slidename=slidename, model=MODEL, simint=SIMINT_DEFAULT
        )


if __name__ == "__main__":
    main()
