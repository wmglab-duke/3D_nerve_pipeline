"""Created on Wed Sep  7 19:58:14 2022.

@author: Daniel
"""

# TODO: make this into a video animation
import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

os.chdir('../..')

for sample, samplenum in zip(["2L", "2R", "3R", "5R", "6L", "6R"], [253, 273, 373, 573, 653, 673]):
    fibers = []
    fiberpath = os.path.join(os.getcwd(), fr'samples\{samplenum}\models\0\sims\3\3D_fiberset')
    slidespath = os.path.join(os.getcwd(), fr'input\slides\{sample}slides.obj')
    # load pickled slidelist
    with open(slidespath, 'rb') as f:
        slidelist = pickle.load(f)

    # for slide in slidelist: #this was a fix for incorrect scaling
    #     slide.scale(0.5)

    # load each fiber file and append to list
    for file in os.listdir(fiberpath):
        if file.endswith('.dat'):
            fibers.append(np.loadtxt(os.path.join(fiberpath, file), skiprows=1))

    # %%     makde video
    os.makedirs('vids', exist_ok=True)

    # importing movie py libraries
    from moviepy.editor import VideoClip
    from moviepy.video.io.bindings import mplfig_to_npimage

    fps = 30

    # duration of the video
    duration = len(slidelist) / fps

    # matplot subplot
    fig, ax = plt.subplots(dpi=200)

    # method to get frames
    def make_frame(t):
        # clear
        ax.clear()

        i = int(t * fps)

        slide = slidelist[i]

        zpos = i * 20  # 20 microns per slice
        # get list of x,y coordinates for each fiber at this z position
        xs = []
        ys = []
        for fiber in fibers:
            # find the index of the closest z value
            idx = (np.abs(fiber[:, 2] - zpos)).argmin()
            # get the x,y coordinates at that index
            x = fiber[idx, 0]
            y = fiber[idx, 1]
            # append to slide list
            xs.append(x)
            ys.append(-y)
        # plot the slide and all fiber points
        ax.scatter(xs, ys, s=3, color='red')
        ax.set_title(f'Slide {i}-zpos{zpos}')
        slide.plot(ax=ax)

        # returning numpy image
        return mplfig_to_npimage(fig)

    # creating animation
    animation = VideoClip(make_frame, duration=duration)

    # saving the clip
    animation.write_videofile(f"vids/{sample}.mp4", fps=fps)
