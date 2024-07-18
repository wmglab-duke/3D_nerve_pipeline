"""Defines plotting functions used for analyzing data.

See ``examples/analysis`` for examples of how to use.
"""

import glob
import json
import os
import pickle
import sys
import warnings
from typing import List, Union

import matplotlib.colorbar as cbar
import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import pandas as pd
import seaborn as sns
from nd_line.nd_line import nd_line
from scipy.interpolate import griddata
from scipy.stats import pearsonr
from shapely.geometry import Point

from src.core import Query, Sample, Simulation
from src.utils import Config, Object


def heatmaps(
    *facetdata,
    data=None,
    ax=None,
    **kwargs,
):
    """Create heatmap for a single axis using the _HeatmapPlotter class.

    To create a single heatmap, call the class directly
    Use a Seaborn ``FacetGrid`` to create multiple heatmaps in one figure, using the ``FacetGrid.map()`` method.
    Note that data cannot be aggregated across n_sims
    (e.g., each call of ``heatmaps()`` must recieve only one threshold per fiber).

    :param facetdata: Recieves data from ``FacetGrid`` if using to plot an array.
    :param data: DataFrame to plot, used if manually passing data.
    :param ax: Axis to plot on.
    :param kwargs: Arguments to be passed to ``_HeatmapPlotter`` class constructor.
    :return: Plotting axis.
    """
    if data is None:
        data = pd.concat(facetdata, axis=1)
    # initialize heatmap plotter and pass in all arguments
    plotter = _HeatmapPlotter(data, **kwargs)
    if ax is None:
        ax = plt.gca()

    plotter.plot(ax)

    return ax


class _HeatmapPlotter:
    """Class used to contruct heatmap plots.

    This class should not be called directly by the user. Rather, the
    user should call the ``heatmaps()`` function, which will pass any
    keyword arguments to this class's constructor.
    """

    def __init__(
        self,
        data,
        mode: str = 'fibers',
        sample_object=None,
        sim_object=None,
        missing_color='red',
        suprathresh_color='blue',
        subthresh_color='green',
        cutoff_thresh=None,
        cmap=None,
        colorbar=True,
        min_max_ticks=False,
        cuff_orientation=False,
        plot_outers=False,
        cbar_kws=None,
        scatter_kws=None,
        line_kws=None,
        min_thresh=None,
        max_thresh=None,
        color=None,
    ):
        """Initialize heatmap plotter.

        :param data: DataFrame containing data to plot.
        :param mode: Plotting mode. There are multiple options:

            * ``'fibers'``: Plot a point for each fiber, using a heatmap of thresholds for color.
            * ``'fibers_on_off'``: Plot a point for each fiber. If the fiber threshold is above cutoff_thresh,
              suprathresh_color is used. Otherwise, subthresh_color is used.
            * ``'inners'``: Plot each inner as filled in, using a heatmap of thresholds for color.
              The mean threshold for that inner is used,
              thus if only one fiber is present per inner, that threshold is used.
            * ``'inners_on_off'``: Plot each inner as filled in. If the mean inner threshold is above cutoff_thresh,
              suprathresh_color is used. Otherwise, subthresh_color is used.

        :param sample_object: Sample object to use for plotting. Automatically loaded if not provided.
        :param sim_object: Simulation object to use for plotting. Automatically loaded if not provided.
        :param missing_color: Color to use for missing data.
        :param suprathresh_color: Color to use for suprathresh data.
        :param subthresh_color: Color to use for subthresh data.
        :param cutoff_thresh: Threshold to use for plotting on_off modes.
        :param cmap: Color map to override default.
        :param colorbar: Whether to add a colorbar.
        :param min_max_ticks: Whether to add only the minimum and maximum ticks to the colorbar.
        :param cuff_orientation: Whether to plot a point for the cuff orientation.
        :param plot_outers: Whether to plot the fascicle outers.
        :param cbar_kws: Keyword arguments to pass to matplotlib.colorbar.Colorbar.
        :param scatter_kws: Keyword arguments to pass to matplotlib.pyplot.scatter.
        :param line_kws: Keyword arguments to pass to matplotlib.pyplot.plot.
        :param min_thresh: Minimum threshold to use for plotting. Use this to override the default minimum.
        :param max_thresh: Maximum threshold to use for plotting. Use this to override the default maximum.
        :param color: Color passed in by seaborn when using FacetGrid. Not used.
        """
        # add variables to self from input args
        self.mappable = None
        self.fiber_colors = self.inner_colors = None
        self.sample_index = self.sim_index = self.model_index = self.n_sim_index = None
        self.plot_outers = plot_outers
        self.cmap = cmap
        self.min_max_ticks = min_max_ticks
        self.colorbar = colorbar
        self.cutoff_thresh = cutoff_thresh
        self.missing_color = missing_color
        self.suprathresh_color = suprathresh_color
        self.subthresh_color = subthresh_color
        self.mode = mode
        self.sample = sample_object
        self.sim = sim_object
        self.color = color
        self.cbar_kws = cbar_kws if cbar_kws is not None else {}
        self.scatter_kws = scatter_kws if scatter_kws is not None else {}
        self.scatter_kws.setdefault('s', 100)
        self.line_kws = line_kws if line_kws is not None else {}
        self.max_thresh = max(data.threshold) if max_thresh is None else max_thresh
        self.min_thresh = min(data.threshold) if min_thresh is None else min_thresh
        self.cuff_orientation = cuff_orientation

        # run setup in preparation for plotting
        self.validate(data)
        self.get_objects()
        self.create_cmap()
        self.determine_colors(data)
        self.data = data

    def plot(self, ax):
        """Make heatmap plot.

        :param ax: Axis to plot on.
        :return: Plotting axis.
        """
        self.set_ax(ax)

        if self.colorbar and self.mode != 'on_off':
            self.add_colorbar(ax)

        if self.cuff_orientation:
            self.plot_cuff_orientation(ax)

        self.plot_inners_fibers(ax)

        return ax

    def plot_inners_fibers(self, ax):
        """Plot inners and fibers using the colors determined in determine_colors().

        :param ax: axis to plot on
        """
        self.sample.slides[0].plot(
            final=False,
            fix_aspect_ratio=True,
            fascicle_colors=self.inner_colors,
            ax=ax,
            outers_flag=self.plot_outers,
            inner_format='k-',
            scalebar=True,
            line_kws=self.line_kws,
        )
        if np.any([bool(x) for x in self.fiber_colors]):
            if self.mode != 'fibermeshgrid':
                self.scatter_kws['c'] = self.fiber_colors
                self.sim.fibersets[0].plot(ax=ax, scatter_kws=self.scatter_kws)
            else:
                points = self.sim.fibersets[0].xy_points()
                # set up meshgrid from x and y points, where values comes from meshgridcolor
                inner_index = 0
                for fascicle in self.sample.slides[0].fascicles:
                    for inner in fascicle.inners:
                        # find indices of points that are in this inner
                        inds = np.where([inner.contains(Point(x, y)) for x, y in points])[0]
                        inner_index += 1
                        # get points to pass in to
                        self.add_inner_meshgrid(inner, np.array(self.fiber_colors)[inds], np.array(points)[inds])

    def add_inner_meshgrid(self, inner, thresholds, points):
        """Add a meshgrid to the plot.

        :param ax: axis to plot on
        :param inner: inner to plot
        :param thresholds: thresholds to use
        :param points: points to use
        """
        # get the points for the inner
        xs = inner.points[:, 0]
        ys = inner.points[:, 1]
        # go through each point in xs and zs and find the closest point in self.data
        # then add the threshold from that index from self.data.threshold to minthreshes
        minthreshes = []
        for i in range(len(xs)):
            minthreshes.append(thresholds[np.argmin(np.sqrt((xs[i] - points[0]) ** 2 + (ys[i] - points[1]) ** 2))])
        minthreshes = np.array(minthreshes)
        # minthreshes will be stepwise, first find the midpoint of each step
        # first find all places where np.diff is nonzero
        diff = np.diff(minthreshes)
        if not np.all(diff == 0):
            diff = np.where(diff != 0)[0]
            # now find the midpoint of each step
            midpoints = []
            for i in range(len(diff)):
                midpoints.append(int((diff[i] + diff[i - 1]) / 2))
            midpoints = np.array(midpoints)
            # replace the original minthreshes with nan where its not a midpoint
            for i in range(len(minthreshes)):
                if i not in midpoints:
                    minthreshes[i] = np.nan
            # now interpolate the nans, connecting the ends of the array
            minthreshes = np.interp(
                np.arange(len(minthreshes)),
                np.where(~np.isnan(minthreshes))[0],
                minthreshes[~np.isnan(minthreshes)],
                period=len(minthreshes),
            )
        x = np.concatenate([points[:, 0], xs])
        y = np.concatenate([points[:, 1], ys])
        z = np.concatenate([thresholds, minthreshes])
        xi = np.linspace(min(x), max(x), 1000)
        yi = np.linspace(min(y), max(y), 1000)
        zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='linear')
        import matplotlib

        clip = matplotlib.path.Path(inner.points[:, :2])
        plt.pcolormesh(
            xi,
            yi,
            zi,
            cmap=self.cmap,
            vmin=self.min_thresh,
            vmax=self.max_thresh,
            clip_path=(clip, plt.gca().transData),
        )  # TODO: figure out how to pass mappable in
        # plt.contour(xi, yi, zi, cmap=self.cmap, levels=15)
        # plt.scatter(x, y, c=z, cmap=self.cmap, vmin=self.min_thresh, vmax=self.max_thresh)

    def create_cmap(self):
        """Create color map and mappable for assigning colorbar and ticks."""
        if self.cmap is None:
            cmap = plt.cm.get_cmap('viridis')
            cmap.set_bad(color='w')
            cmap = cmap.reversed()
            self.cmap = cmap
        mappable = plt.cm.ScalarMappable(
            cmap=self.cmap,
            norm=mplcolors.Normalize(vmin=self.min_thresh, vmax=self.max_thresh),
        )
        self.mappable = mappable

    @staticmethod
    def set_ax(ax):
        """Remove axis elements.

        :param ax: axis to plot on
        """
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')

    def determine_colors(self, threshdf):
        """Determine colors for inners and fibers based on user selected mode.

        :param threshdf: DataFrame of thresholds.
        """

        def _mapthresh(thresh):
            return tuple(self.cmap((thresh - self.min_thresh) / (self.max_thresh - self.min_thresh)))

        inner_count = 0
        for fascicle in self.sample.slides[0].fascicles:
            for _inner in fascicle.inners:
                inner_count += 1

        inner_color_list = []
        fiber_color_list = []
        for inner in range(inner_count):
            # get inner threshold and add the appropriate color to the list
            innerthresh = np.mean(threshdf.query(f'inner=={inner}').threshold)
            if innerthresh is np.nan and self.mode in ['inners', 'inners_on_off']:
                inner_color_list.append(self.missing_color)
                warnings.warn(
                    'Missing at least one fiber threshold, color will appear as missing color (defaults to red).',
                    stacklevel=2,
                )
            elif self.mode == 'inners':
                inner_color_list.append(_mapthresh(innerthresh))
            elif self.mode == 'inners_on_off':
                inner_color_list.append(
                    self.suprathresh_color if innerthresh > self.cutoff_thresh else self.subthresh_color
                )
            else:
                inner_color_list.append(None)
        for fiber_index in pd.unique(threshdf['master_fiber_index']):
            # get fiber threshold and add the appropriate color to the list
            fiberthresh = np.mean(threshdf.query(f'master_fiber_index=={fiber_index}').threshold)
            if fiberthresh is np.nan and self.mode in ['fibers', 'fibers_on_off', 'fibermeshgrid']:
                warnings.warn(
                    'Missing fiber threshold, color will appear as missing color (defaults to red).', stacklevel=2
                )
                fiber_color_list.append(self.missing_color)
            elif self.mode == 'fibers':
                fiber_color_list.append(_mapthresh(fiberthresh))
            elif self.mode == 'fibers_on_off':
                fiber_color_list.append(
                    self.suprathresh_color if fiberthresh > self.cutoff_thresh else self.subthresh_color
                )
            elif self.mode == 'fibermeshgrid':
                fiber_color_list.append(fiberthresh)
            else:
                fiber_color_list.append(None)
        # set colors for inners and fibers
        self.inner_colors, self.fiber_colors = inner_color_list, fiber_color_list

    def add_colorbar(self, ax):
        """Add colorbar to heatmap plot.

        :param ax: axis to plot on
        """
        # set default ticks if not provided
        if 'ticks' not in self.cbar_kws:
            self.cbar_kws['ticks'] = (
                tick.AutoLocator() if not self.min_max_ticks else [self.min_thresh, self.max_thresh]
            )
        # generate colorbar
        cb_label = r'mA'
        cb: cbar.Colorbar = plt.colorbar(mappable=self.mappable, ax=ax, **self.cbar_kws)
        cb.ax.set_title(cb_label)

    def get_objects(self):
        """Get sample and sim objects for plotting."""
        if self.sample is None:
            self.sample = Query.get_object(Object.SAMPLE, [self.sample_index])
        if self.sim is None:
            self.sim = Query.get_object(Object.SIMULATION, [self.sample_index, self.model_index, self.sim_index])

    def validate(self, data):
        """Check that data is valid for plotting.

        :param data: DataFrame of thresholds.
        """
        assert self.mode in ['fibers', 'inners', 'fibers_on_off', 'inners_on_off', 'fibermeshgrid'], 'Invalid mode'
        if self.mode in ['fibers_on_off', 'inners_on_off']:
            assert self.cutoff_thresh is not None, 'Must provide cutoff threshold for on/off mode.'
        # make sure only one sample, model, sim, and nsim for this plot
        for index in ['sample', 'model', 'sim', 'nsim']:
            assert (
                len(pd.unique(data[index])) == 1
            ), f'Only one {index} allowed for this plot. Append something like q.threshold_data.query(\'{index}==0\')'
            setattr(self, index + '_index', pd.unique(data[index])[0])

    def plot_cuff_orientation(self, ax):
        """Plot the orientation of the cuff.

        :param ax: axis to plot on
        """
        # calculate orientation point location (i.e., contact location)
        # get radius of sample
        try:
            r = self.sample.slides[0].nerve.mean_radius()
        except AttributeError:
            r = self.sample.slides[0].fascicles[0].outer.mean_radius()
        # get orientation angle from slide
        if self.sample.slides[0].orientation_angle is None:
            raise ValueError("Cannot plot orientation if orientation angle was not defined in the slide.")
        theta = self.sample.slides[0].orientation_angle
        # load add_ang from model.json cofiguration file
        with open(Query.build_path(Config.MODEL, [self.sample_index, self.model_index])) as f:
            model_config = json.load(f)
        # add any cuff rotation
        theta += np.deg2rad(model_config.get('cuff').get('rotate').get('add_ang'))
        ax.scatter(r * 1.2 * np.cos(theta), r * 1.2 * np.sin(theta), 300, 'red', 'o')


def ap_loctime(  # noqa: C901
    query_object: Query,
    n_sim_filter: List[int] = None,
    plot: bool = False,
    plot_distribution: bool = False,
    n_sim_label_override: str = None,
    model_labels: List[str] = None,
    save: bool = False,
    subplots=False,
    nodes_only=False,
    amp=0,
):
    """Plot time and location of action potential initiation.

    :param query_object: Query object to use for plotting.
    :param n_sim_filter: List of n_sim values to plot.
    :param plot: Whether to plot the ap location node for each fiber.
    :param plot_distribution: Whether to plot action potential initiation node distribution.
    :param n_sim_label_override: Label to use for n_sim.
    :param model_labels: Labels to use for models.
    :param save: Whether to save the plot.
    :param subplots: Whether to plot in subplots.
    :param nodes_only: Whether to plot only nodes.
    :param amp: Amplitude of action potential.
    """
    # loop samples
    for sample_index, sample_results in [(s['index'], s) for s in query_object._result.get('samples')]:
        print(f'sample: {sample_index}')

        # loop models
        for model_index, model_results in [(m['index'], m) for m in sample_results.get('models')]:
            print(f'\tmodel: {model_index}')

            # loop sims
            for sim_index in model_results.get('sims', []):
                print(f'\t\tsim: {sim_index}')

                sim_object = query_object.get_object(Object.SIMULATION, [sample_index, model_index, sim_index])

                if subplots is True:
                    fig, axs = plt.subplots(ncols=len(sim_object.master_product_indices), nrows=2, sharey="row")

                # loop nsims
                for n_sim_index, (potentials_product_index, _waveform_index) in enumerate(
                    sim_object.master_product_indices
                ):
                    print(f'\t\t\tnsim: {n_sim_index}')

                    _, *_, fiberset_index = sim_object.potentials_product[potentials_product_index]

                    # skip if not in existing n_sim filter
                    if n_sim_filter is not None and n_sim_index not in n_sim_filter:
                        print('\t\t\t\t(skip)')
                        continue

                    # directory of data for this (sample, model, sim)
                    sim_dir = query_object.build_path(
                        Object.SIMULATION, [sample_index, model_index, sim_index], just_directory=True
                    )

                    # directory for specific n_sim
                    n_sim_dir = os.path.join(sim_dir, 'n_sims', str(n_sim_index))

                    # directory of fiberset (i.e., points and potentials) associated with this n_sim
                    fiberset_dir = os.path.join(sim_dir, 'fibersets', str(fiberset_index))

                    # the simulation outputs for this n_sim
                    outputs_path = os.path.join(n_sim_dir, 'data', 'outputs')

                    # path of the first inner, first fiber vm(t) data
                    inner = 0
                    if plot_distribution:
                        fiber_indices = len(glob.glob(fiberset_dir + '/*.dat'))
                    else:
                        fiber_indices = [0]
                    ap_nodes = []
                    for fiber_index in fiber_indices:
                        vm_t_path = os.path.join(
                            outputs_path, f'ap_loctime_inner{inner}_fiber{fiber_index}_amp{amp}.dat'
                        )

                        # load vm(t) data (see path above)
                        # each row is a snapshot of the voltages at each node [mV]
                        # the first column is the time [ms]
                        # first row is holds column labels, so this is skipped (time, node0, node1, ...)
                        aploc_data = np.loadtxt(vm_t_path, skiprows=0)

                        aploc_data[np.where(aploc_data == 0)] = float('Inf')

                        time = min(aploc_data)

                        node = np.argmin(aploc_data)
                        ap_nodes.append(node)

                        # create message about AP time and location findings
                        message = f't: {time} ms, node: {node + 1} (of {len(aploc_data) + 2})'
                        if time != float('inf'):
                            print(f'\t\t\t\t\t\t{message}')
                        else:
                            print('No action potential occurred.')
                            continue

                        # plot the AP location with voltage trace
                        # create subplots
                        if plot or save:
                            print(f'\t\t\t\tinner: {inner} \n \t\t\t\t\tfiber: {fiber_index}')
                            if subplots is not True:
                                fig, axes = plt.subplots(1, 1)
                                axes = [axes]
                            else:
                                axes = [axs[0][n_sim_index], axs[1][n_sim_index]]
                            # load fiber coordinates
                            fiber = np.loadtxt(os.path.join(fiberset_dir, f'{fiber_index}.dat'), skiprows=1)
                            nodefiber = fiber[0::11, :]

                            # plot fiber coordinates in 2D
                            if nodes_only is not True:
                                axes[0].plot(fiber[:, 0], fiber[:, 2], 'b.', label='fiber')
                            else:
                                axes[0].plot(nodefiber[:, 0], nodefiber[:, 2], 'b.', label='fiber')

                            # plot AP location
                            axes[0].plot(fiber[11 * node, 0], fiber[11 * node, 2], 'r*', markersize=10)

                            # location display settings
                            n_sim_label = (
                                f'n_sim: {n_sim_index}' if (n_sim_label_override is None) else n_sim_label_override
                            )
                            model_label = '' if (model_labels is None) else f', {model_labels[model_index]}'
                            axes[0].set_xlabel('x location, µm')

                            axes[0].set_title(f'{n_sim_label}{model_label}')
                            if subplots is not True:
                                axes[0].legend(['fiber', f'AP ({message})'])
                            else:
                                axes[0].legend(['fiber', 'AP'])

                            plt.tight_layout()

                            # voltages display settings
                            if subplots is not True or n_sim_index == 0:
                                axes[0].set_ylabel('z location, µm')
                            plt.tight_layout()

                        # display
                        if save:
                            plt.savefig(
                                (
                                    f'out/analysis/ap_time_loc_{sample_index}_{model_index}_{sim_index}_{n_sim_index}_'
                                    f'inner{inner}_fiber{fiber_index}.png'
                                ),
                                dpi=300,
                            )

                        if plot:
                            plt.show()
                    if plot_distribution:
                        total_nodes = len(aploc_data) + 2
                        plt.hist(ap_nodes, range=(0, total_nodes), bins=total_nodes)
                        plt.title('Distribution of AP initiation node sites')
                        plt.xlabel('Node Index')
                        plt.ylabel('Node Count')
                        plt.show()
                        if save:
                            plt.savefig(
                                (f'out/analysis/ap_time_loc_distribution_{sample_index}_{model_index}_{sim_index}.png'),
                                dpi=300,
                            )


def _get_object(mode: Object, indices: List[int]) -> Union[Sample, Simulation]:
    """Load pickled python object from file.

    :param mode: Mode of object to load.
    :param indices: Indices of object to load.
    :return: Object of specified mode and indices.
    """
    with open(_build_path(mode, indices), 'rb') as obj:
        return pickle.load(obj)


def _build_path(
    mode: Union[Config, Object],
    indices: List[int] = None,
) -> str:
    """Build path to pickled python object.

    :param mode: Mode of object to load.
    :param indices: Indices of object to load.
    :return: Path to pickled python object.
    """
    result = ''

    if indices is None:
        indices = [
            0,
            0,
            0,
        ]  # dummy values... will be stripped from path later bc just_directory is set to True

    if mode == Config.SAMPLE:
        result = os.path.join('samples', str(indices[0]), 'sample.json')
    elif mode == Config.MODEL:
        result = os.path.join('samples', str(indices[0]), 'models', str(indices[1]), 'model.json')
    elif mode == Config.SIM:
        result = os.path.join('config', 'user', 'sims', f'{indices[0]}.json')
    elif mode == Object.SAMPLE:
        result = os.path.join('samples', str(indices[0]), 'sample.obj')
    elif mode == Object.SIMULATION:
        result = os.path.join(
            'samples',
            str(indices[0]),
            'models',
            str(indices[1]),
            'sims',
            str(indices[2]),
            'sim.obj',
        )

    return result


def datamatch_agg(dest, dat3d, importval, merge=False):
    dest[importval + '3d'] = np.nan  # todo: update this to take input columns to use for matching
    for i in range(len(dest)):
        row = dest.iloc[i, :]
        val = dat3d[
            (dat3d["model"] == row['model'])
            & (dat3d["sim"] == row['sim'])
            & (dat3d["nerve_label"] == row['nerve_label'])
            & (dat3d["nsim"] == row['nsim'])
            & (dat3d["level"] == row['level'])
        ][importval]
        val = list(val)
        if len(val) != 1:
            sys.exit('issue here')
        dest.iloc[i, -1] = val[0]
    if np.any(dest[importval] == np.nan):
        sys.exit('issue here too')
    return dest


def datamatch_merge(dest, dat3d, importval, merge_cols):
    merged = dest.merge(dat3d[merge_cols + [importval]], on=merge_cols, how="left", suffixes=["", "3d"])

    if merged[importval + "3d"].isna().any():
        raise RuntimeError('Found nan value after merge.')

    return merged


def datamatch(dest, dat3d, importval, merge=False):
    dest[importval + '3d'] = np.nan
    for i in range(len(dest)):
        row = dest.iloc[i, :]
        val = dat3d[
            (dat3d["model"] == row['model'])
            & (dat3d["sim"] == row['sim'])
            & (dat3d["nerve_label"] == row['nerve_label'])
            & (dat3d["nsim"] == row['nsim'])
            & (dat3d["master_fiber_index"] == row['master_fiber_index'])
        ][importval]
        val = list(val)
        if len(val) == 0:
            raise RuntimeError('No match for master fiber index.')
        elif len(val) > 1:
            raise RuntimeError('Found more than one match for master fiber index.')
        dest.iloc[i, -1] = val[0]
    if np.any(dest[importval] == np.nan):
        sys.exit('issue here too')
    return dest


def datamatchlist(dest, dat3d, importvals, merge=False):
    for ival in importvals:
        dest[ival + '3d'] = np.nan
        for i in range(len(dest)):
            row = dest.iloc[i, :]
            val = dat3d[
                (dat3d["model"] == row['model'])
                & (dat3d["sim"] == row['sim'])
                & (dat3d["nerve_label"] == row['nerve_label'])
                & (dat3d["nsim"] == row['nsim'])
                & (dat3d["master_fiber_index"] == row['master_fiber_index'])
            ][ival]
            val = list(val)
            if len(val) != 1:
                sys.exit('issue here')
            dest.iloc[i, -1] = val[0]
        if np.any(dest[ival] == np.nan):
            sys.exit('issue here too')
    return dest


def rename_var(df, di):
    for variable, values in di.items():
        for old, new in values.items():
            df = df.replace(to_replace={variable: old}, value=new)
    return df


def plot_colorthresh(samp2d, samp3d, model, simdex, nerve_label):
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp2d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat2d = q.data()
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat3d = q.data(source_sample=samp2d)
    dat2d = datamatch(dat2d, dat3d, 'threshold')
    # %% Renaming
    redict = {
        # "nsim": {
        #     0: 'fiber diameter: 2\u03BCm',
        #     1: 'fiber diameter: 5\u03BCm',
        #     2: 'fiber diameter: 8\u03BCm',
        #     3: 'fiber diameter: 11\u03BCm',
        #     4: 'fiber diameter: 13\u03BCm',
        # }
    }
    dat2d = dat2d.rename(columns={'threshold': '2D', 'threshold3d': '3D'})
    datre = rename_var(dat2d, redict)
    dat2dnew = datre.drop(columns='3D').rename(columns={'2D': 'threshold'})
    dat2dnew['dataset'] = '2D'
    dat3dnew = datre.drop(columns='2D').rename(columns={'3D': 'threshold'})
    dat3dnew['dataset'] = '3D'
    datfinal = pd.concat([dat2dnew, dat3dnew], sort=True)
    # datre = dat2d
    # %%
    sns.set(font_scale=1.5)
    plotdata = datfinal[datfinal['sample'] == samp2d]
    g = sns.catplot(
        data=plotdata,
        kind='swarm',
        col='nsim',
        hue='inner',
        y='threshold',
        x='dataset',
        sharey=False,
        palette='colorblind',
    )
    plt.subplots_adjust(top=0.85)
    plt.suptitle(f'Activation thresholds by fascicle (Sample {nerve_label}, 2D slice: {samp2d})')
    axs = g.axes.ravel()
    axs[0].set_ylabel('Activation threshold (mA)')
    plt.subplots_adjust(top=0.85)
    for i, ax in enumerate(g.axes.ravel()):
        corr = {}
        thisdat = dat2d[(dat2d["nsim"] == i) & (dat2d["sample"] == samp2d)]
        corr[samp2d] = round(pearsonr(thisdat['2D'], thisdat['3D'])[0], 3)
        leg = ax.legend(
            labels=["r=" + str(corr[samp2d])], handlelength=0, handletextpad=0, fancybox=True, loc='lower center'
        )
        for item in leg.legendHandles:
            item.set_visible(False)
    plt.savefig(f'out/analysis/{simdex}/colorthresh{nerve_label}-{samp2d}.png', dpi=500)


def plot_colorthreshstrip(samp2d, samp3d, model, simdex, nerve_label):
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp2d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat2d = q.data()
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat3d = q.data(source_sample=samp2d)
    dat2d = datamatch(dat2d, dat3d, 'threshold')
    # %% Renaming
    redict = {
        # "nsim": {
        #     0: 'fiber diameter: 2\u03BCm',
        #     1: 'fiber diameter: 5\u03BCm',
        #     2: 'fiber diameter: 8\u03BCm',
        #     3: 'fiber diameter: 11\u03BCm',
        #     4: 'fiber diameter: 13\u03BCm',
        # }
    }
    dat2d = dat2d.rename(columns={'threshold': '2D', 'threshold3d': '3D'})
    datre = rename_var(dat2d, redict)
    dat2dnew = datre.drop(columns='3D').rename(columns={'2D': 'threshold'})
    dat2dnew['dataset'] = '2D'
    dat3dnew = datre.drop(columns='2D').rename(columns={'3D': 'threshold'})
    dat3dnew['dataset'] = '3D'
    datfinal = pd.concat([dat2dnew, dat3dnew], sort=True)
    # datre = dat2d
    # %%
    sns.set(font_scale=1.5)
    plotdata = datfinal[datfinal['sample'] == samp2d]
    g = sns.catplot(
        data=plotdata,
        kind='strip',
        col='nsim',
        hue='inner',
        y='threshold',
        x='dataset',
        sharey=False,
        dodge=True,
        palette='colorblind',
    )
    plt.subplots_adjust(top=0.85)
    plt.suptitle(f'Activation thresholds by fascicle (Sample {nerve_label}, 2D slice: {samp2d})')
    axs = g.axes.ravel()
    axs[0].set_ylabel('Activation threshold (mA)')
    plt.subplots_adjust(top=0.85)
    for i, ax in enumerate(g.axes.ravel()):
        corr = {}
        thisdat = dat2d[(dat2d["nsim"] == i) & (dat2d["sample"] == samp2d)]
        corr[samp2d] = round(pearsonr(thisdat['2D'], thisdat['3D'])[0], 3)
        leg = ax.legend(
            labels=["r=" + str(corr[samp2d])], handlelength=0, handletextpad=0, fancybox=True, loc='lower center'
        )
        for item in leg.legendHandles:
            item.set_visible(False)
    plt.savefig(f'out/analysis/{simdex}/colorthresh{nerve_label}-{samp2d}.png', dpi=500)


def get_peri_site(samp3d, samp2d, model, simdex, nerve_label, source_sim=3):
    """Go through each fiber threshold and find the perineurium thickness at the activation site."""
    q3 = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [0], 'sim': [simdex]},
        }
    ).run()
    dat3d = q3.data(source_sample=samp2d)
    dat3d['threed'] = True
    dat3z = get_actual_zpos(dat3d, samp3d, model, simdex, source_sim=source_sim, xy=True)
    dat3z['peri_thk_act_site'] = np.nan
    dat3z.rename(columns={'threshold': 'threshold3d'}, inplace=True)
    # Loop through each fiber and find the perineurium thickness at the activation site
    for i in range(0, len(dat3z)):
        # get the z position of the activation site
        zpos = dat3z.loc[i, 'activation_zpos']
        # Load in pickled slidelist
        with open(f'input/slides/{nerve_label}slides.obj', 'rb') as f:
            slidelist = pickle.load(f)
        # find the slice nearest the activation site.
        slice_spacing = 20  # microns
        slice_index = int(round(zpos / slice_spacing))
        # Get the slide
        slide = slidelist[slice_index]
        slide.scale(1.2)  # shrinkage correction
        slide.scale(0.5)  # wrong scaling correction #TODO remove this
        # Use the x and y activation pos to create a shapely point, then find which inner from
        # the slide contains that point
        point = Point(dat3z.loc[i, 'activation_xpos'], dat3z.loc[i, 'activation_ypos'])
        inner = None
        try:
            inner = [inner for fasc in slide.fascicles for inner in fasc.inners if inner.contains(point)][0]
        except IndexError:
            print('ope')
            iteration = 0
            innersave = [inner for fasc in slide.fascicles for inner in fasc.inners]
            innerlist = [x.deepcopy() for x in innersave]
            while inner is None:
                if iteration > 5:
                    plt.figure()
                    plt.scatter(dat3z.loc[i, 'activation_xpos'], dat3z.loc[i, 'activation_ypos'])
                    [inner.plot() for inner in innerlist]
                    # slide.plot()
                    plt.show()
                    plt.title(
                        f'slide_index: {slice_index}\nzpos-{zpos}\nmaster_fiber{dat3z.loc[i, "master_fiber_index"]}'
                    )
                    break
                else:
                    iteration += 1
                [x.offset(distance=10) for x in innerlist]
                try:
                    inner = innersave[int(np.where([inner.contains(point) for inner in innerlist])[0])]
                    print('ope fixed')
                except IndexError:
                    pass
        if inner is not None:
            fit = {'a': 0.03702, 'b': 10.5}
            thk = fit.get("a") * 2 * np.sqrt(inner.area() / np.pi) + fit.get("b")
            dat3z.loc[i, 'peri_thk_act_site'] = thk
    return dat3z


def get_datamatch(samples2d, samp3d, model, simdex, nerve_label, tortuosity=False, source_sim=None, cuffspan=None):
    global ax
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': samples2d, 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat2d = q.data(tortuosity=tortuosity)
    q3 = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat3d = q3.data(source_sample=samples2d[0], tortuosity=tortuosity, source_sim=source_sim, cuffspan=cuffspan)
    dat2d = datamatch(dat2d, dat3d, 'threshold')
    if tortuosity:
        dat2d = datamatch(dat2d, dat3d, 'tortuosity')
    dat2d['nerve_label'] = nerve_label
    return dat2d


def corrcalc(data, comparison):
    corrs = []
    for nsim in pd.unique(data['nsim']):
        # ax.set_title(f'fiber diam: {s}μm')
        corr = {}
        for sample in pd.unique(data['sample']):
            thisdat = data[(data["nsim"] == nsim) & (data["sample"] == sample)]
            corr[sample] = round(pearsonr(thisdat[comparison[0]], thisdat[comparison[1]])[0], 3)
        corrs.append(corr)
    return corrs


def plot_correlation(samples2d, samp3d, model, simdex, nerve_label):
    global ax
    corrs = []
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': samples2d, 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat2d = q.data()
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat3d = q.data(source_sample=samples2d[0])
    dat2d = datamatch(dat2d, dat3d, 'threshold')
    # %%
    import seaborn as sns
    from scipy.stats import pearsonr

    sns.set_theme()
    sns.set(font_scale=1.5)
    dat2d = dat2d.rename(columns={'sample': 'Slice'})
    # %%
    g = sns.lmplot(
        data=dat2d,
        x="threshold",
        y="threshold3d",
        hue="Slice",
        height=5,
        col='nsim',
        facet_kws={'sharey': False, 'sharex': False},
    )
    axs = g.axes.ravel()
    axs[0].set_ylabel('3D threshold (mA)')
    plt.suptitle(f'Activation threshold correlation for sample {nerve_label}', fontsize=25)
    plt.subplots_adjust(top=0.85, right=0.93)
    new_labels = ['Anodic\nLeading', 'Cathodic\nLeading']
    for t, l in zip(g._legend.texts, new_labels):
        t.set_text(l)
    for i, ax in enumerate(g.axes.ravel()):
        # ax.set_title(f'fiber diam: {s}μm')
        corr = {}
        for sample in samples2d:
            thisdat = dat2d[(dat2d["nsim"] == i) & (dat2d["Slice"] == sample)]
            corr[sample] = round(pearsonr(thisdat['threshold'], thisdat['threshold3d'])[0], 3)
        ax.legend(labels=["r=" + str(corr[sample]) for sample in samples2d])
        ax.set_xlabel('2D threshold (mA)')
        corrs.append(corr)
    g.savefig(f'out/analysis/{simdex}/threshcorr_{nerve_label}', dpi=400)
    return corrs


def plot_colorjoint(samp2d, samp3d, model, simdex, nerve_label):
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp2d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat2d = q.data()
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat3d = q.data(source_sample=samp2d)
    dat2d = datamatch(dat2d, dat3d, 'threshold')
    # %%
    sns.set(font_scale=1.5)
    dat2d = dat2d[dat2d['sample'] == samp2d]
    for nsim in pd.unique(dat2d.nsim):
        plotdata = dat2d.query(f'nsim=={nsim}')
        g = sns.jointplot(data=plotdata, x="threshold3d", y="threshold", hue="inner", palette='colorblind')
        ax = g.ax_joint
        idat = plotdata
        min_thresh = min([min(idat.threshold), min(idat.threshold3d)])
        max_thresh = max([max(idat.threshold), max(idat.threshold3d)])
        limits = (min_thresh, max_thresh)
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.plot(limits, limits, color='red')
        plt.savefig(f'out/analysis/{simdex}/colorjoint{nerve_label}-{samp2d}-{nsim}.png', dpi=500)


def plot_oneone(samples2d, samp3d, model, simdex, nerve_label):
    global ax
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': samples2d, 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat2d = q.data()
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat3d = q.data(source_sample=samples2d[0])
    dat2d = datamatch(dat2d, dat3d, 'threshold')
    # %%
    sns.set_theme()
    sns.set(font_scale=1.5)
    dat2d = dat2d.rename(columns={'sample': 'Slice'})
    # %%
    g = sns.lmplot(
        data=dat2d,
        x="threshold3d",
        y="threshold",
        hue="Slice",
        height=5,
        col='nsim',
        facet_kws={'sharey': False, 'sharex': False},
    )
    axs = g.axes.ravel()
    axs[0].set_ylabel('2D threshold (mA)')
    plt.suptitle(f'Activation threshold correlation for sample {nerve_label}', fontsize=25)
    plt.subplots_adjust(top=0.85, right=0.93)
    new_labels = ['Anodic\nLeading', 'Cathodic\nLeading']
    for t, l in zip(g._legend.texts, new_labels):
        t.set_text(l)
    for i, ax in enumerate(g.axes.ravel()):
        # ax.set_title(f'fiber diam: {s}μm')
        idat = dat2d.query(f'nsim=={i}')
        min_thresh = min([min(idat.threshold), min(idat.threshold3d)])
        max_thresh = max([max(idat.threshold), max(idat.threshold3d)])
        limits = (min_thresh, max_thresh)
        corr = {}
        for sample in samples2d:
            thisdat = dat2d[(dat2d["nsim"] == i) & (dat2d["Slice"] == sample)]
            corr[sample] = round(pearsonr(thisdat['threshold'], thisdat['threshold3d'])[0], 3)
        ax.legend(labels=["r=" + str(corr[sample]) for sample in samples2d])
        ax.set_xlabel('3D threshold (mA)')
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.plot(limits, limits, color='red')
    g.savefig(f'out/analysis/{simdex}/oneone_{nerve_label}', dpi=400)


def plot_dose_response(samples2d, samp3d, model, simdex, nerve_label):
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': samples2d, 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat2d = q.data()
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat3d = q.data(source_sample=samples2d[0])
    data = pd.concat([dat2d, dat3d])
    data.reset_index(inplace=True)
    for i in range(len(pd.unique(data['nsim']))):
        # ax = axs[i]
        plt.figure()
        sns.color_palette("tab10")
        # plt.figure()
        plotdata = data[data.nsim == i]
        plotdata = plotdata[plotdata['sample'] != 670]
        sns.ecdfplot(data=plotdata, x='threshold', hue='sample', palette='colorblind')
        plt.xscale('log')
        plt.ylabel('Proportion of Fibers Activated')
        plt.xlabel('Activation Threshold (mA, log scale)')
        plt.title(f'Threshold eCDF for {nerve_label}')
        plt.savefig(f'out/analysis/{simdex}/{nerve_label}_{i}_ecdf.png', dpi=400, bbox_inches='tight')


def ap_plot(samp2d, samp3d, model, simdex, cuff_contacts, source_sim=None):
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp2d], 'model': [0], 'sim': [simdex]},
        }
    ).run()
    dat2d = q.data()
    dat2d['threed'] = False
    q3 = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [0], 'sim': [simdex]},
        }
    ).run()
    dat3d = q3.data(source_sample=samp2d)
    dat3d['threed'] = True
    # %%
    dat3z = get_actual_zpos(dat3d, samp3d, model, simdex, source_sim=source_sim)
    dat2d = dat2d.rename(columns={'long_ap_pos': 'activation_zpos'})
    apdat = pd.concat([dat3z, dat2d], sort=True)
    # %%
    redict = {"sample": {samp2d: '2D', samp3d: '3D'}}
    datre = rename_var(apdat, redict)
    datre.loc[:, 'activation_zpos'] = datre.loc[:, 'activation_zpos'] / 1000
    # %%
    g = sns.catplot(
        x="sample",
        y='activation_zpos',
        hue='sample',
        col="nsim",
        data=datre,
        kind='strip',
        height=5,
        aspect=0.4,
        linewidth=0,
        order=['2D', '3D'],
        sharey=True,
    )
    axs = g.axes
    axs[0][0].set_ylabel('Activation z-position (mm)')
    # for i, s in enumerate([2, 5, 8, 11, 13]):
    #     axs[0][i].set_title(f'fiber diam: {s}μm')
    plt.subplots_adjust(top=0.88)
    plt.suptitle(f"Activation z-position for sample {samp2d}", fontsize=15)
    for ax in axs[0]:
        ln = ax.hlines([cuff_contacts], 0.4, 0.6, color='red')
    axs[0][-1].legend([ln], ['Cuff Contacts'], loc='center right')
    g.savefig(f'out/analysis/{simdex}/{samp2d}_zpos.png', dpi=400)


def get_actual_zpos(dat3d, samp3d, model, sim, source_sim=None, xy=False):
    fiberdir = os.path.join(
        'samples',
        str(samp3d),
        'models',
        str(model),
        'sims',
        str(sim if source_sim is None else source_sim),
        '3D_fiberset',
    )
    fibers3d = [x for x in os.listdir(fiberdir) if x.endswith('.dat')]
    for file in fibers3d:
        f_ind = int(os.path.splitext(file)[0])
        coord = np.loadtxt(os.path.join(fiberdir, file), skiprows=1)
        fiberline = nd_line(coord)
        datnew = dat3d[(dat3d["model"] == model) & (dat3d["sim"] == sim) & (dat3d["master_fiber_index"] == f_ind)]
        for index, row in datnew.iterrows():
            actual_zpos = fiberline.interp(row['long_ap_pos'])[2]
            dat3d.loc[index, 'activation_zpos'] = actual_zpos
            if xy:
                dat3d.loc[index, 'activation_xpos'] = fiberline.interp(row['long_ap_pos'])[0]
                dat3d.loc[index, 'activation_ypos'] = fiberline.interp(row['long_ap_pos'])[1]
    return dat3d
