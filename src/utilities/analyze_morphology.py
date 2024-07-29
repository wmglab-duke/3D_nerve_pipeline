import glob
import json
import math
import os
import pickle
import re
import sys

import cv2
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from shapely.ops import unary_union

sys.path.append('../../..')
sys.path.append('../../../Scripts')
from src.core import Fascicle, Nerve, Slide, Trace
from src.utils import ContourMode, DeformationMode, DownSampleMode, MaskSpaceMode, NerveMode
from threedclass import FascicleConnectivityMap

root = os.path.abspath(os.path.join(__file__, '..', '..', '..'))

ascentdir = r'D:\threed_final\input\slides'

os.makedirs(ascentdir, exist_ok=True)
mpl.rcParams['figure.dpi'] = 400


def get_sorted_image_list(path, pattern):
    # uses a regex patter to get a sorted list of images (sorted by the number in the patter)
    x = glob.glob(path + '/' + pattern)
    x = [os.path.split(el)[-1] for el in x]
    x = [os.path.splitext(el)[0] for el in x]
    x = [re.split(r'\D', el) for el in x]
    x = [[el for el in lis if el != ''] for lis in x]
    x = [el[-1] for el in x]
    x = sorted(x, key=int)
    check = [int(el) for el in x]
    start = min(check)
    stop = max(check) + 1
    if check != list(range(start, stop)):
        sys.exit("error")
    x = [pattern.replace('*', el) for el in x]
    return x, start, stop


def get_slide_dims(slidelist):
    allbd = [slide.bounds() for slide in slidelist]
    allbd = np.array(allbd)
    allbd = (min(allbd[:, 0]), min(allbd[:, 1]), max(allbd[:, 2]), max(allbd[:, 3]))
    dims = (
        (math.floor(allbd[0]) - buffer, math.ceil(allbd[2]) + buffer),
        (math.floor(allbd[1]) - buffer, math.ceil(allbd[3]) + buffer),
    )
    dims = [[abs(y) for y in x] for x in dims]
    dims = tuple([(-max(x), max(x)) for x in dims])
    return dims


detailed = True

pickling = False

movie = False

pickleonly = False

nerve_mode = NerveMode.PRESENT
deform_mode = DeformationMode.PHYSICS
os.chdir(os.path.join(root))

# samples = ['2LDS5','2RDS5','3RDS5','5RDS5','6LDS5','6RDS5','2LDS5def', '3RDS5def', '5RDS5def', '6RDS5def','2LDS5_fine','2LDS5_im','2LDS5_up','2LDS5_down']
# samples = ['2LDS5']
samples = ['2LDS5', '2RDS5', '3RDS5', '5RDS5', '6LDS5', '6RDS5']

if detailed:
    fig, axs = plt.subplots(nrows=5, ncols=len(samples), sharex=True)

fmaps = {}

morphdir = os.path.join(root, 'morphology')
os.makedirs(morphdir, exist_ok=True)
fasccounts = {}
diamcounts = {}

for sample_index, sample in enumerate(samples):
    config3d_path = os.path.join(root, 'config', f'{sample}.json')

    with open(config3d_path) as f:
        params = json.load(f)

    dire = os.path.join(root, 'datanew', sample, params['path']['slides'])

    dire = dire.replace('/work/wmglab/dpm42/3d_vns/datanew/', 'D:/work/threed/datanew/')

    buffer = 0

    n, start, stop = get_sorted_image_list(dire + '/i/', '*.tif')
    imgs = get_sorted_image_list(dire + '/i/', '*.tif')[0]

    print('Generating mapslides...')

    slides = []

    for i in range(start, stop):
        if i % 100 == 0:
            print(100 * i / stop)
        img_nerve = cv2.imread(dire + f'/n/{n[i]}', -1)

        if len(img_nerve.shape) > 2 and img_nerve.shape[2] > 1:
            img_nerve = img_nerve[:, :, 0]
        contour, _ = cv2.findContours(np.flipud(img_nerve), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        nerve = Nerve(Trace([point + [0] for point in contour[0][:, 0, :]]))

        fascicles = []

        fascicles = Fascicle.to_list(dire + f'/i/{imgs[i]}', None, ContourMode.NONE, MaskSpaceMode.IMAGE)

        slide: Slide = Slide(
            fascicles,
            nerve,
            nerve_mode,
            will_reposition=(deform_mode != DeformationMode.NONE),
        )

        # shift slide about (0,0)
        slide.move_center(np.array([0, 0]), target=params["preprocess"]["center_target"])
        slide.scale(5)  # TODO: make this not hard-coded, this is the um per pixel, also to be correct should be 5
        # slide.smooth_traces(params["preprocess"]["nsmoothing"], params["preprocess"]["fsmoothing"])
        slide.move_center(np.array([0, 0]), target=params["preprocess"]["center_target"])
        slide.generate_perineurium(fit={'a': 0.03702, 'b': 10.50})
        slide.z = i * 20  # TODO make this not hard coded, slice spacing
        slides.append(slide)
    fasccounts[f'{sample}'] = [len(s.fascicles) for s in slides]
    zsep = 20
    span = [1300, 6700]
    centerspan = span - np.mean(span)
    midpoint = int(stop - start / 2)
    # using zsep as microns per slice, find the slice matching the top and bottom of centerspan (which is in microns)
    zspan = [int(np.round(centerspan[0] / zsep)), int(np.round(centerspan[1] / zsep))]
    finalspan = [midpoint + zspan[0], midpoint + zspan[1]]
    allpoly = unary_union([s.nerve.polygon() for s in slides])
    points = list(allpoly.exterior.coords)
    alltrace = Trace(points)
    alltrace.down_sample(DownSampleMode.KEEP, int(np.floor(alltrace.points.size / 100)))
    nerve_copy = alltrace
    x, y, r_bound = nerve_copy.make_circle()
    print(f'Bounding radius for {sample} is {r_bound}')
    if pickling:
        with open(os.path.join(ascentdir, sample + 'slides.obj'), 'wb') as f:
            pickle.dump(slides, f)
    if movie:
        # importing movie py libraries
        from moviepy.editor import VideoClip
        from moviepy.video.io.bindings import mplfig_to_npimage

        dims = get_slide_dims(slides)
        # duration of the video
        duration = 50

        # matplot subplot
        fig, ax = plt.subplots(dpi=300)

        # method to get frames
        def make_frame(i):
            # clear
            ax.clear()

            # plotting line
            ax.set_ylim([dims[1][0], dims[1][1]])
            ax.set_xlim([dims[0][0], dims[0][1]])
            slides[int(i * 50)].plot(ax=ax)

            # returning numpy image
            return mplfig_to_npimage(fig)

        # creating animation
        animation = VideoClip(make_frame, duration=duration)

        # displaying animation with auto play and looping
        animation.write_videofile(os.path.join(morphdir, f'{sample}vid.mp4'), fps=50)

    if pickleonly:
        continue

    fmaps[sample] = FascicleConnectivityMap(slides)
    fmaps[sample].generate_map()
    # %%
    if detailed:
        areanum, nervediamnum, fascdiamnum, fasccountnum, mergesplitnum = range(5)

        sns.set(font_scale=1, style='ticks', context='paper')
        fig, axs = plt.subplots(1, 5, sharey=True)
        plt.subplots_adjust(wspace=0.2)
        i_a = []
        n_a = []
        i_n = []
        i_as = []
        p_ia = []
        i_alla = []

        for slide in slides:
            n_a.append(slide.nerve.area())
            areas = [i.area() for f in slide.fascicles for i in f.inners]
            i_a.append(sum(areas))
            i_n.append(len(areas))
            i_as.append(np.mean(areas))
            i_alla.append(areas)
            p_ia.append(sum(areas) / slide.nerve.area())

        xpos = (np.array(range(len(slides))) * 20) / 1000

        branches = [np.array(fmaps[sample].merge_inds) * 20 / 1000, np.array(fmaps[sample].split_inds) * 20 / 1000]

        xs = []
        yds = []
        for i in range(len(slides))[::100]:
            ys = [inner.ecd() for f in slides[i].fascicles for inner in f.inners]
            yds.extend(ys)
            for y in ys:
                xs.append(i * 20 / 1000)

        diamcounts[f'{sample}'] = []
        for i in range(len(slides))[::1]:
            ys = [inner.ecd() for f in slides[i].fascicles for inner in f.inners]
            diamcounts[f'{sample}'].extend(ys)

        axs[0].set_ylabel('Position along nerve (cm)')
        axs[fascdiamnum].scatter(np.array(yds) / 1000, xs, marker='+', color='k')
        # axs[0].set_xscale('log')
        axs[fascdiamnum].set_xlabel('Fascicle diameter\n(mm)')
        axs[areanum].plot(np.array(n_a) * 1e-6, xpos, color='k', label='Whole Nerve')
        axs[fasccountnum].plot(i_n, xpos, color='k')
        axs[nervediamnum].set_xlabel('Nerve diameter\n(mm)')
        axs[nervediamnum].plot(np.sqrt(np.array(n_a) / np.pi) * 2e-3, xpos, color='k')
        axs[fasccountnum].set_xlabel('Fascicle count')
        axs[areanum].plot(np.array(i_a) * 1e-6, xpos, color='k', ls=':', label='Endoneurium')
        axs[areanum].set_xlabel(r'Area ($\mathrm{mm^2}$)')
        # axs[areanum].legend()
        axs[mergesplitnum].eventplot(branches, orientation='vertical', color=['k'])
        fig.set_size_inches(14, 6)
        # print(ax.get_xlim())
        axs[mergesplitnum].set_xticks([0, 1])
        axs[mergesplitnum].set_xlim([-1, 2])
        axs[mergesplitnum].set_xticklabels(["merge", "split"], rotation=0)
        # axs[4].legend(labels=['merges','splits'],bbox_to_anchor=[2.5,1])
        # plt.suptitle(f"Morphological Data for {sample[:2]}")
        # plt.subplots_adjust(hspace=0, wspace=0.3)
        fig.set_size_inches(8, 3, forward=True)
        fig.subplots_adjust(top=0.9)
        # fig.savefig(os.path.join(morphdir, 'detailed.png'), dpi=400, bbox_inches='tight')
        plt.ylim(reversed(plt.ylim()))
        plt.yticks(np.arange(0, 60, 10), [int(x) for x in reversed(np.arange(0, 60, 10) / 10)])
        sns.despine()
        for ax in axs[:-1]:
            ax.set_xlim([0, None])
            difftick = np.diff(ax.get_xticks())[0]
            modtick = np.diff(ax.get_xlim())[0] % difftick
            ax.set_xlim([None, ax.get_xlim()[1] + difftick - modtick])

        # Add text with arrows
        axs[areanum].annotate(
            'Endoneurium',
            xy=(0.15, 1),
            xytext=(0.15, 1.3),
            xycoords='axes fraction',
            textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', lw=1.5, color='k'),
            ha='center',
            va='center',
        )
        axs[areanum].annotate(
            'Whole\nnerve',
            xy=(0.45, 1),
            xytext=(0.45, 1.15),
            xycoords='axes fraction',
            textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', lw=1.5, color='k'),
            ha='center',
            va='center',
        )
    # sys.exit()
# %% new plots

import seaborn as sns

sns.set(font_scale=1, style='white', context='paper')
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 400
bc = [fmap.branch_count for fmap in fmaps.values()]
fc = [len(fmap.map_data) for fmap in fmaps.values()]
mf = [np.mean([len(s.fascicles) for s in fmap.slides]) for fmap in fmaps.values()]

plt.figure()
plt.scatter(mf, bc, color='k')
plt.xlabel('Mean fascicle count')
plt.ylabel('Branch count')
# plt.title('Branch count vs mean fascicle count')
# add best fit line and print correlation coefficient
m, b = np.polyfit(mf, bc, 1)
plt.plot(mf, m * np.array(mf) + b, color='black', label=f'$R^2$={np.corrcoef(mf,bc)[0,1]**2:.2f}')
plt.legend()
plt.gcf().set_size_inches(3, 2)
plt.savefig(os.path.join(morphdir, 'bcounts.png'), dpi=400)
# %%
# sns.reset_orig()
plt.figure()


fc = pd.DataFrame(fasccounts)
maxcount = np.amax(np.amax(fc))
# remove DS5 from each column name
fc.columns = [col.replace('DS5', '') for col in fc.columns]

g = sns.histplot(fc, element='poly', bins=range(maxcount + 1), palette='colorblind', fill=False, legend=False)
plt.xlim([0, 20])
plt.xlabel('Number of fascicles')
plt.ylabel('Count')
plt.xticks([0, 5, 10, 15, 20])
# remove legend title
# hand,lab = g.get_legend_handles_labels()
# sns.move_legend(g,[1,0.3],title='Nerve',ncol=1,frameon=False)
plt.gcf().set_size_inches(3, 2)
sns.despine()
# change dataframe to long form for seaborn
fcnew = fc.stack().reset_index()
fcnew.columns = ['slice', 'sample', 'count']
# plot with y axis as slice, x axis as count, and color as sample
plt.figure()
sns.lineplot(data=fcnew, y='count', x='slice', hue='sample', palette='colorblind')
plt.ylabel('Fascicle Count')
plt.xlabel('Slice')
plt.yticks([0, 10, 20])
plt.gca().get_legend().set_title('')
plt.gcf().set_size_inches(4, 3)
plt.legend(ncol=3)

# redo using matplotlib with x axis as slice, y axis as count, and color as sample
plt.figure()
# get seaborn colorblind palette
palette = sns.color_palette('colorblind')
for sample in fc.columns:
    plt.plot(fc[sample], np.array(fc.index) * 20e-3, label=sample, color=palette[fc.columns.get_loc(sample)])
plt.ylabel('Position (mm)')
plt.xlabel('Fascicle Count')
plt.legend()
plt.gcf().set_size_inches(3, 4)
plt.xticks([0, 10, 20])

# %%
# diamcounts is a dict of lists of sizes foreach sample. Each sample has different number of sizes
# convert to long form for seaborn (so each row is a size and a sample)
alldata = []
for sample in diamcounts.keys():
    for size in diamcounts[sample]:
        alldata.append([sample, size])
fs = pd.DataFrame(alldata)
fs.columns = ['nerve', 'sizes']
fs.nerve = fs.nerve.str.replace('DS5', '')
plt.figure()
g = sns.histplot(
    fs,
    x='sizes',
    hue='nerve',
    element='poly',
    palette='colorblind',
    fill=False,
    legend=True,
    bins=np.arange(0, 1600, 100),
)
sns.move_legend(g, [1.1, 0.1], title='')
plt.xlabel('Fascicle Diameter (μm)')
plt.gcf().set_size_inches(3, 2)
sns.despine()
