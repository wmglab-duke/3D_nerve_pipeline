"""Takes segmented nerve/fascicle, preprocesses the images,
creates the perineurium, generates a connectivity map, and creates fibersets."""

# edit to use threedmodel TODO max importance
# also may need to move over cuff models from 3D? to see #TODO
import argparse
import gc
import glob
import json
import math
import os
import pickle
import re
import shutil
import sys
import time
import traceback
from copy import deepcopy

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import splev, splprep
from scipy.ndimage import binary_fill_holes
from shapely.ops import unary_union
from skimage import morphology

scriptloc = os.path.abspath(__file__)
root = os.path.split(os.path.split(scriptloc)[0])[0]
scriptroot = os.path.split(scriptloc)[0]
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "1"
sys.path.insert(0, os.path.abspath(root))
sys.path.insert(0, os.path.join(root, 'src', 'utilities'))
from nd_line.nd_line import nd_line
from src.core import Deformable, Fascicle, FiberSet, Model, Nerve, Slide, Trace
from src.runner import Runner
from src.utils import (
    Config,
    Configurable,
    ContourMode,
    DeformationMode,
    DownSampleMode,
    MaskSpaceMode,
    NerveMode,
    ReshapeNerveMode,
    SetupMode,
)
from threedclass import FascicleConnectivityMap

parser = argparse.ArgumentParser()
parser.add_argument("config_script", nargs="?")
parser.add_argument("-q", "--quit-premesh", action="store_true")
parser.add_argument("-A", "--run-all-steps", action="store_true")
parser.add_argument("-p", "--mesher-progress", action="store_true")
args = parser.parse_args()

if args.config_script is not None:
    configname = args.config_script
else:
    configname = '5RDS5def_MCT.json'
config3d_path = os.path.join('..', 'config', configname)
# %% defs
pseudonym = os.path.splitext(os.path.split(configname)[-1])[0]
print(f'Running {pseudonym}')


def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
    return


def run_scanip(configpath, configdict, script, progress=False):
    with open(configpath, 'w') as f:
        json.dump(configdict, f, indent=2)
    os.chdir(env['SIMPLEWARE_SCANIP_PATH'])
    os.chdir(env['SIMPLEWARE_SCANIP_PATH'])
    executable = "./ConsoleScanIP" if sys.platform.startswith("linux") else "consolescanip.exe"
    progstring = "--no-progress" if not progress else ''
    os.system(
        f'{executable} {progstring} --force-queuing'
        f' --disable-undo --run-script="{root}/scripts/{script}.py" --input-value="{configpath}"'
    )
    os.chdir(os.path.split(os.path.split(scriptloc)[0])[0])
    print('returned to python')


def get_sorted_image_list(imgpath, pattern):
    # uses a regex pattern to get a sorted list of images (sorted by the number in the patter)
    x = glob.glob(imgpath + '/' + pattern)
    x = [os.path.split(el)[-1] for el in x]
    x = [os.path.splitext(el)[0] for el in x]
    x = [re.split(r'\D', el) for el in x]
    x = [[el for el in lis if el != ''] for lis in x]
    x = [el[-1] for el in x]
    x = sorted(x, key=int)
    check = [int(el) for el in x]
    if len(check) == 0:
        raise FileNotFoundError("No images found in " + imgpath)
    begin = min(check)
    end = max(check) + 1
    if check != list(range(begin, end)):
        raise ValueError("Some images are missing")
    x = [pattern.replace('*', el) for el in x]
    if len(x) == 0:
        raise FileNotFoundError("No images found")
    return x, begin, end


class AscentConfigs(Configurable):
    def __init__(self):
        Configurable.__init__(self)


with open(config3d_path) as f:
    params = json.load(f)

modelfile = params['config']['model']

samplefile = 'samples' + '/' + params['config']['sample']

ascent_configs = AscentConfigs()

ascent_configs.add(SetupMode.NEW, Config.MODEL, root + f'/config/models/{modelfile}')
ascent_configs.add(SetupMode.NEW, Config.SIM, root + '/config/sim.json')
ascent_configs.add(SetupMode.NEW, Config.SAMPLE, root + f'/config/{samplefile}')
ascent_configs.add(SetupMode.NEW, Config.FIBER_Z, root + '/config/fiber_z.json')

deformer = ascent_configs.configs[Config.SAMPLE.value]['modes']['deform']

if deformer == "PHYSICS":
    deform_mode = DeformationMode.PHYSICS
elif deformer == "NONE":
    deform_mode = DeformationMode.NONE
else:
    sys.exit('invalid deform mode')

with open(root + '/config/system/env.json') as f:
    env = json.load(f)

nerve_mode = NerveMode.PRESENT

projpath = os.path.join(env['ASCENT_PROJECT_PATH'], 'datanew', os.path.split(config3d_path)[-1].replace('.json', ''))

if 'project_path' in params:
    print(
        f'''Project path will be overwritten with new smart pathing...
        From {params['project_path']}
        to {projpath}
        '''
    )
params['project_path'] = projpath

with open(config3d_path, 'w') as f:
    json.dump(params, f, indent=2)

for key, value in params['path'].items():
    params['path'][key] = params['project_path'] + '/' + value

buffer = params['preprocess']['buffer']

# TODO break point config

preproc = False

imfills = False

slidegen = False

rundeform = False

geometry = False

mesh = False

model = False

fibergen = False

extract = True

skip_getpots = True

quit_premesh = False

no_remesh = True

skipsave = True

# TODO: make each model source from datanew directory instead of main model directory

quit_premesh = args.quit_premesh or quit_premesh

if args.run_all_steps:
    preproc = imfills = slidegen = geometry = mesh = model = fibergen = extract = True
    skipsave = False
    if args.quit_premesh:
        quit_premesh = True
    else:
        quit_premesh = False

if no_remesh and os.path.exists(params['path']['mesh'] + '/mesh.nas'):
    preproc = imfills = slidegen = geometry = mesh = False

resample_factor = params['preprocess']['upsample_z']

for path in params['path'].values():
    ensure_dir(path)

if not os.path.exists(params['project_path'] + '/model.json'):
    shutil.copyfile(root + f'/config/models/{modelfile}', params['project_path'] + '/model.json')
# %% Image preprocessing
os.chdir(scriptroot)
sipsource = params['sourcedir']
infile = params['sourcefile']
sourcesipfile = os.path.join(sipsource, infile)
if not os.path.exists(sourcesipfile):
    raise FileNotFoundError(f'source .sip file not found: {sourcesipfile}')

# start of "3D" run
print('Preprocessing images...')

os.makedirs(os.path.join(sipsource, "preprocessed"), exist_ok=True)
preprocpath = sipsource + '/preprocessed/' + os.path.splitext(infile)[0]
# skip preprocessing if the modified date of the preprocessed folder is newer than the sip file
if os.path.exists(preprocpath) and os.path.getmtime(preprocpath) > os.path.getmtime(sipsource + '/' + infile):
    print('preprocessed images found, skipping preprocessing')
    preproc = False

if preproc:
    # set up scanIP config and run sip file
    sipconfig = {
        'outpath': os.path.abspath(preprocpath),
        'infile': os.path.abspath(sipsource + '/' + infile),
        'resample_factor': resample_factor,
        'resample_xy': params['input']['um_per_px'] / params['output']['image_um_per_px'],
        'um_per_px': params['input']['um_per_px'],
        'um_per_slice': params['input']['um_per_slice'],
        'n_mask': params['input']['n_mask'],
        'i_mask': params['input']['i_mask'],
    }
    scriptname = 'sip_pp_destep'
    sip_config_path = root + '/config/system/configdict.json'
    run_scanip(sip_config_path, sipconfig, scriptname)
# edit um per slice for resampling
params['input']['um_per_slice'] = params['input']['um_per_slice'] / resample_factor
params['input']['um_per_px'] = params['output']['image_um_per_px']

# %% Generate slides
os.chdir(root)
os.chdir('scripts')

n_imgs_pp, _, _ = get_sorted_image_list(preprocpath + '/n', '*.tiff')

i_imgs_pp, start, stop = get_sorted_image_list(preprocpath + '/i', '*.tiff')

if not preproc:
    imin = cv2.imread(preprocpath + '/i/' + i_imgs_pp[0])
print('Filling holes in the images...')


def fill_img_holes(filldir):
    for file in os.listdir(filldir):
        if file.endswith('.tiff') or file.endswith('.tif'):
            img = cv2.imread(os.path.join(filldir, file), cv2.IMREAD_GRAYSCALE)
            img = binary_fill_holes(img)
            img = morphology.remove_small_objects(img, 20)
            imgout = (255 * (img / np.amax(img))).astype(np.uint8)
            cv2.imwrite(os.path.join(filldir, file), imgout)


if preproc or imfills:
    for imtype in ['i', 'n']:
        path = preprocpath + f'/{imtype}'
        fill_img_holes(path)

# make a sample from ASCENT
print('Generating slides...')


def slidegenerator(nerveimg, fascicleimg, scale='all'):
    img_nerve = cv2.imread(nerveimg, -1)
    if len(img_nerve.shape) > 2 and img_nerve.shape[2] > 1:
        img_nerve = img_nerve[:, :, 0]
    contour, _ = cv2.findContours(np.flipud(img_nerve), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    nerve = Nerve(Trace([point + [0] for point in contour[0][:, 0, :]]))
    fascicles = Fascicle.to_list(fascicleimg, None, ContourMode.NONE, MaskSpaceMode.CARTESIAN)
    slide: Slide = Slide(
        fascicles,
        nerve,
        nerve_mode,
        will_reposition=True,
    )
    # shift slide about (0,0)
    slide.move_center(np.array([0, 0]))
    slide.scale(params["input"]["um_per_px"])
    if scale == 'all':
        # shrinkage correction
        slide.scale(1.2)  # TODO this should come from config
    slide.smooth_traces(params["preprocess"]["nsmoothing"], params["preprocess"]["fsmoothing"])
    slide.move_center(np.array([0, 0]))
    if scale == 'all':
        slide.generate_perineurium(fit={'a': 0.03702, 'b': 10.50})
    slide.validate(intersection_target='inners')
    return slide


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


if slidegen:
    slides = []
    # TODO make from ASCENT
    for i in range(start, stop):
        try:
            slide = slidegenerator(preprocpath + f'/n/{n_imgs_pp[i]}', preprocpath + f'/i/{i_imgs_pp[i]}')
            slides.append(slide)
        except Exception:
            slide.plot()
            plt.show()
            traceback.print_exc()
            sys.exit(f'Error generating slide index {i}')
    dims = get_slide_dims(slides)
    nervesave = deepcopy([s.nerve for s in slides])
# %%
print('Generating fascicle connectivity map...')
if slidegen and not skipsave:
    fmap = FascicleConnectivityMap(slides)
    # newreadded
    fmap.generate_map()
    # endnewreadded
    print('writing images...')
    fmap.save_images(
        params['path']['slides'] + '/multi_cleaned',
        dims,
        separate=False,
        print_ids=True,
        resize_factor=1 / params["output"]["image_um_per_px"],
        id_font=params["output"]["id_font"],
    )
    fmap.nervestats(preprocpath + '/nerve_ecd_stats.json')
    # newstatblock
    branchstats = {'fascicule count': len(fmap.map_data), 'branch count': fmap.branch_count}
    savepath = params['path']['slides'] + '/branchstats.json'
    with open(savepath, 'w') as f:
        json.dump(branchstats, f, indent=2)

# %% Generate comsol file


def get_override_r(slidelist):
    return np.mean([np.sqrt(s.nerve.area() / np.pi) for s in slidelist])


def compute_medium_params(configmodel: dict, nslice: int, um_per_slice: float):
    for med in configmodel["medium"]:
        configmodel["medium"][med]["length"] = nslice * um_per_slice + 2 * configmodel["medium"][med]["buffer"]
        configmodel["medium"][med]["shift"] = {'z': -configmodel["medium"][med]["buffer"], 'x': 0, 'y': 0}
    return configmodel


def run_comsol(threed_config_path, runtype):
    print('to java')
    os.chdir(root)
    runner = Runner(0)
    runner.add(SetupMode.OLD, Config.ENV, env)
    runner.populate_env_vars()
    core_name = 'ModelWrapper'

    runner.handoff(threed_config_path, run_type=runtype, class_name=core_name, modelfolder="threedmodel")
    print('returned to python')


def get_minbound_r(slidelist):
    return np.amax([s.nerve.make_circle()[-1] for s in slidelist])


# TODO remove this class
class AscentUtil(Configurable):
    def __init__(self, Configurable):
        # initialize Configurable super class
        Configurable.__init__(self)

    # TODO make from ASCENT
    def compute_cuff_shift_3D(self, modelcon: dict, slidesec: list, rparam: str):  # noqa: N802
        deform_ratio = self.search(Config.SAMPLE, 'deform_ratio')

        if len(slidesec) == 0:
            sys.exit("no slides?")

        if deform_ratio != 1:
            allpoly = unary_union([s.nerve.polygon() for s in slidesec])
            points = list(allpoly.exterior.coords)
            alltrace = Trace(points)
            alltrace.down_sample(DownSampleMode.KEEP, int(np.floor(alltrace.points.size / 100)))
            nerve_copy = alltrace
        else:
            max_index = np.argmax([slide.nerve.area() for slide in slidesec])
            override = np.mean([slide.nerve.area() for slide in slidesec])
            override_r = (override / np.pi) ** 0.5
            nerve_copy = slidesec[max_index].nerve.to_circle(override_r=override_r)
        defslide: Slide = Slide(
            [],
            nerve_copy,
            nerve_mode,
            will_reposition=(deform_mode != DeformationMode.NONE),
        )
        os.chdir('..')
        model = Model()
        model.add(SetupMode.NEW, Config.MODEL, os.path.join(params['project_path'], 'model.json'))
        model.add(SetupMode.NEW, Config.SIM, root + '/config/sim.json')
        with open(root + f'/config/{samplefile}') as f:
            sample_config = json.load(f)
        sample_config['Morphology'] = {}
        sample_config['Morphology']['Nerve'] = {}
        sample_config['Morphology']['Nerve']['area'] = nerve_copy.area()
        model.compute_cuff_shift(
            defslide, sample_config, addl_cuff_buffer=params["output"]["addl_cuff_buffer"], ignore_uncentered=True
        )
        new_model_config = model.configs[Config.MODEL.value]
        modelcon[rparam] = new_model_config['min_radius_enclosing_circle']
        os.chdir('scripts')
        return modelcon

    def deform_sep(self, slidelist, override_sep=None):
        # deforms all slides in slidelist, if mean_r is true will deformo all to the mean radius
        # find a way to check for too much self intersection of the underlying geometry
        self.add(SetupMode.NEW, Config.SAMPLE, params['config']['path'] + f'/{samplefile}')

        for slide in slidelist:
            morph_count = 10

            deform_ratio = 0

            sep_fascicles = override_sep or self.search(Config.SAMPLE, "boundary_separation", "fascicles")
            sep_nerve = None

            if 'nerve' in self.search(Config.SAMPLE, 'boundary_separation'):
                sep_nerve = self.search(Config.SAMPLE, 'boundary_separation', 'nerve')
            sep_nerve = sep_nerve - sep_fascicles / 2

            storenerve = slide.nerve.deepcopy()
            slide.nerve.offset(distance=-sep_nerve)
            slide.nerve.scale(1)

            deformable = Deformable.from_slide(slide, ReshapeNerveMode.CIRCLE)

            movements, rotations = deformable.deform(
                morph_count=morph_count,
                render=False,
                minimum_distance=sep_fascicles,
                ratio=deform_ratio,
                progress_bar=False,
            )

            for move, angle, fascicle in zip(movements, rotations, slide.fascicles):
                fascicle.shift(list(move) + [0])
                fascicle.rotate(angle)

            slide.nerve = storenerve
            # shift slide about (0,0)
            slide.move_center(np.array([0, 0]))
        return slidelist

    # TODO make it so that you have to ok image and also make text bigger and white border and also pixel target should be the center of the text
    def deform(
        self,
        runstring,
        span,
        target_radius,
        transition_distance,
        outpath,
        startrad=5000,
        epipad=25,
    ):
        # set up scanIP config and run sip file
        sipconfig = {
            'span': span,
            'epi_radius': target_radius,
            'epipad': epipad,
            'infile': preprocpath + '/debug_prebinary.sip',
            'out': outpath,
            'transition_distance': transition_distance,
            'zcenter': 0.5 * (params["input"]["um_per_slice"] * (span[1] - span[0]) + 2 * transition_distance) / 1000,
            'cuffspan': params["input"]["um_per_slice"] * (span[1] - span[0]) / 1000,  # todo: these are wrong now
            'start_radius': startrad / 1000,
            'outpath': outpath,
        }
        sip_config_path = outpath + '/sipdefconfig.json'
        scriptname = 'sip_deformer_p'
        if not os.path.exists(outpath + '/predeform.stl'):
            run_scanip(sip_config_path, sipconfig, scriptname)
        if not os.path.exists(outpath + '/postdeform.stl'):
            run_comsol(config3d_path, runstring)
        if not os.path.exists(outpath + '/postdeform.stl'):
            sys.exit("postdeform.stl was not created, please manually generate and rerun")
        if not os.path.exists(outpath + '/postdeformed.sip'):
            scriptname = 'sip_postdeform'
            run_scanip(sip_config_path, sipconfig, scriptname)
        return sipconfig

    def deforms_linear(self, slidelist, dest_r, stopspot=None, end_deform_ratio=None, rlist=None, spline=True):
        # linearly deforms along slidelist, moving from deform_ratio=1 and dest_r to deform_ratio=0 and the slide r
        self.add(SetupMode.NEW, Config.SAMPLE, params['config']['path'] + f'/{samplefile}')

        morph_count = 100

        master_deform_ratio = self.search(Config.SAMPLE, 'deform_ratio')
        sep_fascicles = self.search(Config.SAMPLE, "boundary_separation", "fascicles")
        sep_nerve = self.search(Config.SAMPLE, "boundary_separation", "nerve") - sep_fascicles / 2

        for i, slide in enumerate(slidelist):
            ratio = (len(slidelist) - i) / (len(slidelist))
            if end_deform_ratio is not None:
                deform_ratio = (master_deform_ratio - end_deform_ratio) * ratio + end_deform_ratio
            else:
                deform_ratio = master_deform_ratio * ratio

            pre_area = slide.nerve.area()
            if rlist is not None:
                assert len(slidelist) == len(rlist)
                final_r = rlist[i]
            else:
                final_r = np.sqrt(pre_area / np.pi) * (1 - deform_ratio) + dest_r * deform_ratio

            slide.nerve.offset(distance=-sep_nerve)
            slide.nerve.scale(1)
            slide.nerve.points = np.flip(slide.nerve.points, axis=0)  # set points to opencv orientation
            if spline:
                points = np.array(slide.nerve.points)

                tck, u = splprep(points.T, u=None, s=0.0, per=1, k=1)
                u_new = np.linspace(u.min(), u.max(), 250)
                x_new, y_new, z0 = splev(u_new, tck, der=0)
                slide.nerve = Trace(np.stack([x_new, y_new, z0], axis=1))
                slide.nerve.shift([-slide.nerve.centroid()[0], -slide.nerve.centroid()[1], 0])
                assert np.all(np.isclose(slide.nerve.centroid(), [0, 0]))

            deformable = Deformable.from_slide(slide, ReshapeNerveMode.CIRCLE, override_r=final_r - sep_nerve)

            partially_deformed_nerve = Deformable.deform_steps(
                deformable.start, deformable.end, morph_count, deform_ratio
            )[-1]

            if deform_ratio != 1 and partially_deformed_nerve is not None:
                partially_deformed_nerve.shift(-np.asarray(list(partially_deformed_nerve.centroid()) + [0]))
                slide.nerve = partially_deformed_nerve
                slide.nerve.offset(distance=sep_nerve)
            else:
                slide.nerve = slide.reshaped_nerve(ReshapeNerveMode.CIRCLE, override_r=final_r)
            new_area = np.pi * final_r**2
            slide.nerve.scale((new_area / slide.nerve.area()) ** 0.5)

            # shift slide about (0,0)
            slide.move_center(np.array([0, 0]))
            assert np.all(np.isclose(slide.nerve.centroid(), [0, 0]))

        return slidelist


# initialize util
os.chdir(root)
util = AscentUtil(Configurable)
util.add(SetupMode.NEW, Config.MODEL, os.path.join(params['project_path'], 'model.json'))
util.add(SetupMode.NEW, Config.SIM, root + '/config/sim.json')
util.add(SetupMode.NEW, Config.SAMPLE, root + f'/config/{samplefile}')

os.chdir('scripts')
# %%Apply deformation if specified.
if slidegen:
    with open(params['project_path'] + '/model.json') as f:
        model_config = json.load(f)
    # if 'LivaNova' not in model_config['cuff']['preset']:
    # print('WARNING: only livanova modified to work with 3D cuff export, results on other cuffs may vary')
    # deform and get r for zspan1
    transition_distance = 2500  # um #TODO: from config
    transition_n = math.ceil(transition_distance / params['input']['um_per_slice'])
    cuff_config: dict = util.load(os.path.join(root, "config", "system", "cuffs", model_config['cuff']['preset']))
    zcuff = model_config['cuff']['shift']['z']
    # get zspan and calculate slide section
    zspan1 = np.array(cuff_config['zspan1'])
    midpoint = zcuff / params["input"]["um_per_slice"] + (len(slides) + 1) / 2
    slidespan1 = [int(x / params["input"]["um_per_slice"] + midpoint) for x in zspan1]
    slide_section1 = slides[math.floor(slidespan1[0]) : math.ceil(slidespan1[1])]
    r1 = get_override_r(slide_section1)
    deformspan1 = (
        int(slidespan1[0] - transition_distance / params["input"]["um_per_slice"]),
        int(slidespan1[1] + transition_distance / params["input"]["um_per_slice"]),
    )
    deform_section1 = slides[deformspan1[0] : deformspan1[1]]
    if any(num < 0 or num > len(slides) - 1 for num in slidespan1):
        raise RuntimeError("slide out of bounds")
    model_config = util.compute_cuff_shift_3D(model_config, slide_section1, 'min_radius_enclosing_circle1')

    # get zspan and calculate slide section
    if cuff_config.get('zspan2') is not None:
        zspan2 = np.array(cuff_config['zspan2'])
        slidespan2 = [int(x / params["input"]["um_per_slice"] + midpoint) for x in zspan2]
        if any(num < 0 or num > len(slides) - 1 for num in slidespan2):
            raise RuntimeError("slide out of bounds")
        slide_section2 = slides[math.floor(slidespan2[0]) : math.ceil(slidespan2[1])]
        r2 = get_override_r(slide_section2)
        deformspan2 = (
            int(slidespan2[0] - transition_distance / params["input"]["um_per_slice"]),
            int(slidespan2[1] + transition_distance / params["input"]["um_per_slice"]),
        )
        deform_section2 = slides[deformspan2[0] : deformspan2[1]]
        model_config = util.compute_cuff_shift_3D(model_config, slide_section2, 'min_radius_enclosing_circle2')
    with open(params['project_path'] + '/model.json', 'w') as f:
        json.dump(model_config, f, indent=2)
    # %%
    # Run deformation algorithm
    if deform_mode == DeformationMode.PHYSICS and not skipsave and rundeform:
        if False:
            assert 'Livanova' in model_config['cuff']['preset'], "Must use Livanova cuff for physics-based deformation"
        if not os.path.exists(params['path']['slides'] + '/doubledeform/deformed.sip'):
            print('newdefcode')
            newslides = util.deform_sep(deepcopy(slides), override_sep=30)
            print("Generating output images...")
            dimthis = get_slide_dims(newslides)
            thisfmap = FascicleConnectivityMap(newslides)

        def calc_anglediff(target, reference):
            deformedangle = target.ellipse()[-1]
            ref_angle = reference.ellipse()[-1]
            anglediff = deformedangle - ref_angle
            if np.abs(anglediff) > 90:
                print('diff', anglediff)
                anglediff = np.sign(anglediff) * (np.abs(anglediff) - 180)
                print('def', deformedangle)
                print('ref', ref_angle)
                print('diffnow', anglediff)
            return anglediff

        if cuff_config.get('zspan2') is not None:
            print("running double deform algorithm. Do NOT use for any cuff other than LivaNova")
            if not os.path.exists(params['path']['slides'] + '/doubledeform/deformed.sip'):
                thisfmap.save_images(
                    params['path']['slides'] + '/doubledeform',
                    dimthis,
                    separate=True,
                    print_ids=False,
                    resize_factor=1 / params["output"]["image_um_per_px"],
                    id_font=params["output"]["id_font"],
                )

            sipconfigdef = util.deform(
                "doubledeform",
                [slidespan1[0], slidespan2[1]],
                [r1, r2],  # todo fixup target radius (java side change)
                transition_distance,
                params['path']['slides'] + '/doubledeform',
                startrad=get_minbound_r(slide_section1)
                + 1000,  # TODO figure out why this isnt workin (shouldnt need to add 1000)
            )
        else:
            raise NotImplementedError("Currently deformation only supports two cuff spans")
    print('Returned to main script')

    # %%
    def watershed_slide(
        de, active_fascicles, ref, plot=False, blind=False, compactness=0, lastref=None, saveloc=None, interactive=False
    ):  # TODO: use offset instead of shrink for watershed
        def prep_points(points):
            # adjusts plot points to dimensions and formats for PIL
            points = (points - dim_min + buffer)[:, 0:2].astype(int)
            points = tuple(zip(points[:, 0], points[:, 1]))
            return points

        def plot_watershed(interactive=False):
            fig, axes = plt.subplots(ncols=3, figsize=(9, 3), sharex=True, sharey=True)
            ax = axes.ravel()

            ax[0].imshow(np.flipud(image), cmap=plt.cm.gray)
            ax[0].set_title('Overlapping objects')
            if not interactive:
                ax[0].scatter(coords[:, 1], ax[0].get_ylim()[0] - coords[:, 0], s=20)
            # otherwise,print the label value instead of a plot
            else:
                for coord in coords:  # TODO get this working with the flipped image coordinate
                    ax[0].text(
                        coord[1], ax[0].get_ylim()[0] - coord[0], labels[-coord[0], coord[1]], color='k', fontsize=6
                    )
            ax[1].imshow(np.flipud(-distance), cmap=plt.cm.gray)
            ax[1].set_title('Distances')
            ax[2].imshow(labels, cmap=plt.cm.nipy_spectral)
            ax[2].set_title('Separated objects')

            fig.tight_layout()
            if saveloc is not None:
                fig.savefig(saveloc)
                plt.close()
            if plot:
                plt.show()

        # convert fascicles to binary image and use watershed to segment
        from PIL import Image, ImageDraw
        from scipy import ndimage as ndi
        from skimage.feature import peak_local_max
        from skimage.segmentation import watershed

        buffer = 0
        dim_min = [min(x) for x in dims]
        dim = [max(x) for x in dims]
        imdim = [dim[0] + abs(dim_min[0]) + buffer * 2, dim[1] + abs(dim_min[1]) + buffer * 2]
        imgi = Image.new('1', imdim)
        draw = ImageDraw.Draw(imgi)
        for fascicle in de.fascicles:
            if len(fascicle.inners) != 0:
                for inner in fascicle.inners:
                    draw.polygon(prep_points(inner.points[:, 0:2]), fill=1)
            else:
                draw.polygon(prep_points(fascicle.outer.points[:, 0:2]), fill=1)
        image = np.array(imgi).astype(np.uint8) * 255
        distance = ndi.distance_transform_edt(image)
        if blind:
            coords = peak_local_max(distance, min_distance=50, labels=image)
        # get coords from active fascicles
        else:
            coords = np.array(
                prep_points(np.array([f['position'] for f in active_fascicles.values()]).astype(np.int64))
            )
            coords = np.flip(coords, axis=1)
        if lastref is not None:  # use last slide as reference for watershed
            imgi = Image.new('1', imdim)
            draw = ImageDraw.Draw(imgi)
            for fascicle in lastref.fascicles:
                for inner in fascicle.inners:
                    innnn = inner.deepcopy()
                    innnn.scale(0.5)
                    draw.polygon(prep_points(innnn.points[:, 0:2]), fill=1)
            labelimage = np.array(imgi).astype(np.uint8) * 255
            markers, _ = ndi.label(labelimage)
        else:
            mask = np.zeros(distance.shape, dtype=bool)
            mask[tuple(coords.T)] = True
            markers, _ = ndi.label(mask)
        labels = np.flipud(watershed(-distance, markers, mask=image, compactness=compactness))
        if plot or saveloc is not None:
            if not interactive:
                plot_watershed()
            else:
                while True:
                    plot_watershed(interactive=True)
                    get = ''
                    if len(np.unique(labels)) - 1 != len(ref.fascicles):
                        get = input("Choose labels to combine (integers separated by spaces).")
                    if get == '':
                        break
                    # get list of labels to combine
                    labels_to_combine = [int(x) for x in get.split(" ")]
                    # combine labels
                    for label in labels_to_combine:
                        labels[labels == label] = labels_to_combine[0]

        imgn = Image.new('1', imdim)
        draw = ImageDraw.Draw(imgn)
        draw.polygon(prep_points(de.nerve.points[:, 0:2]), fill=1)
        img_nerve = np.flipud(np.array(imgn).astype(np.uint8) * 255)
        if len(img_nerve.shape) > 2 and img_nerve.shape[2] > 1:
            img_nerve = img_nerve[:, :, 0]
        contour, _ = cv2.findContours(np.flipud(img_nerve), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        nerve = Nerve(Trace([point + [0] for point in contour[0][:, 0, :]]))

        fascicles = []
        for label in np.unique(labels):
            if label == 0:
                continue
            contourimg = ((labels == label).astype(np.uint8) * 255).astype(np.uint8)
            contours, _ = cv2.findContours(np.flipud(contourimg), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
            traces = [Trace([item + [0] for item in contour[:, 0, :]]) for contour in contours]
            if len(traces) != 1:
                print("using largest trace, be careful")
                traces = [traces[np.argmax([trace.area() for trace in traces])]]
            outer = traces[0]
            fascicle = Fascicle(outer=outer)
            fascicles.append(fascicle)

        slide: Slide = Slide(
            fascicles,
            nerve,
            nerve_mode,
            will_reposition=(deform_mode != DeformationMode.NONE),
        )
        # shift slide about (0,0)
        slide.move_center(np.array([0, 0]))
        try:
            slide.scale(params["input"]["um_per_px"])
        except Exception as e:
            plt.title('Issue with watershed')
            slide.plot()
            plot_watershed()
            raise e
        # shrinkage correction
        slide.scale(1.2)  # TODO this should come from config
        slide.smooth_traces(params["preprocess"]["nsmoothing"], params["preprocess"]["fsmoothing"])
        slide.move_center(np.array([0, 0]))
        de = slide  # todo: add these steps to making the deformed slide
        de.scale(np.sqrt(ref.nerve.area() / de.nerve.area()))
        return de

    import math

    def rotate(origin, point, angle):
        """Rotate a point counterclockwise by a given angle around a given origin.

        The angle should be given in radians.
        """
        ox, oy = origin
        px, py = point

        qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
        qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
        return qx, qy

    from copy import deepcopy

    slidesave = deepcopy(slides)

    # %%
    def deform_postfix_area_only(deform_section, ref_section, plot=False):
        assert len(deform_section) == len(ref_section)
        defmap = FascicleConnectivityMap(ref_section)
        defmap.generate_map()
        refmap = defmap.get_simplemap()
        begin = [x for x, v in refmap.items() if v['startz'] == 0]
        active_fascicles = {b: {'position': defmap.map_data[b]['points'][0][:2]} for b in begin}
        for i in range(len(deform_section)):
            gc.collect()
            print(i)
            # check for merge or split
            branch = np.any([value['endz'] == i - 1 for value in refmap.values()])
            if branch:  # TODO: make image processing open if split close if merge
                if plot:
                    plt.figure()
                # find the fascicle that is splitting or merging
                splitmerge = [key for key, value in refmap.items() if value['endz'] == i - 1]
                # find where the new fascicle location(s) should be
                plt.title(i)
                deform_section[i - 1].plot(line_kws={"c": "red"}, outers_flag=False, final=False)
                for fascicle in splitmerge:
                    tos = defmap.map_data[fascicle]['to']
                    print(tos)
                    # need to used the deformed rotation to inform this
                    for to in tos:
                        new_position = (
                            np.array(defmap.map_data[to]['points'][0][:2])
                            - np.array(defmap.map_data[fascicle]['points'][-1][:2])
                            + np.array(active_fascicles[fascicle]['position'])
                        )
                        # plot arrow from old to new
                        if plot:
                            plt.arrow(
                                active_fascicles[fascicle]['position'][0],
                                active_fascicles[fascicle]['position'][1],
                                new_position[0] - active_fascicles[fascicle]['position'][0],
                                new_position[1] - active_fascicles[fascicle]['position'][1],
                                width=30,
                            )
                        # need to change this with respect to the fascicle rotation.
                        # so rotate around fascicle centroid with the amount of rotation
                        # being the diff between the ellipses.
                        active_fascicles[to] = {'position': tuple(new_position)}
                    # remove the now gone fascicle
                    active_fascicles.pop(fascicle)
                if plot:
                    plt.title(i)
                    deform_section[i - 1].plot(line_kws={"c": "red"}, outers_flag=False)
                    plt.show()
            # get the slide for this iteration
            ref, de = ref_section[i], deform_section[i]
            desave = deepcopy(de)
            try:
                assert len(ref.fascicles) == len(de.fascicles)
            except AssertionError:
                test = watershed_slide(de, active_fascicles, ref, blind=True, plot=True)
                assert len(ref.fascicles) == len(test.fascicles)
                de = test
                deform_section[i] = de
            # check that all active fascicles are within a fascicle in the deformed slide
            nearests = {}
            neardists = {}
            for num, fdat in active_fascicles.items():
                # find nearest fascicle in deformed
                fasciclepos = fdat['position']
                deformedfascicles = [np.array([f.centroid()[0], f.centroid()[1]]) for f in de.fascicles]
                deformedfascicles = np.array(deformedfascicles)
                dists = np.linalg.norm(deformedfascicles - fasciclepos, axis=1)
                minidx = np.argmin(dists)
                nearests[num] = minidx
                neardists[num] = dists[minidx]
            for num, fdat in active_fascicles.items():  # this code is wrong
                # find closest fascicle to position in reference slide
                original_fascicle_position = defmap.map_data[num]['points'][i - refmap[num]['startz']][:2]
                # get fascicle object matching position from reference slide
                centroids = [np.array([f.centroid()[0], f.centroid()[1]]) for f in ref.fascicles]
                centroids = np.array(centroids)
                dists = np.linalg.norm(centroids - original_fascicle_position, axis=1)
                minidx = np.argmin(dists)
                fdat['position'] = de.fascicles[nearests[num]].centroid()[:2]
                deffasc = de.fascicles[nearests[num]].outer
                reffasc = ref.fascicles[minidx].inners[0]
                deffasc.scale(reffasc.ecd() / deffasc.ecd())
            de.generate_perineurium(fit={'a': 0.03702, 'b': 10.50})
            if plot:
                plt.figure()
                plt.title(f'Slide {i}\nred=corrected size\ngreen=deformed unfixed\nblue=original segmentation')
                ref.plot(final=False, outers_flag=False, line_kws={'c': 'blue', 'label': 'original'})
                de.plot(final=False, outers_flag=False, line_kws={'c': 'red', 'label': 'corrected'})
                desave.plot(final=False, outers_flag=False, line_kws={'c': 'green', 'label': 'deformed'})
                plt.show()
            assert np.all(np.isclose(de.nerve.centroid(), [0, 0]))
        return deform_section

    # %%
    # update new fascicle positions to slides
    if deform_mode == DeformationMode.PHYSICS and rundeform:
        slides = deepcopy(slidesave)
        print('Deform postfix')
        defsection = slides[deformspan1[0] : deformspan2[1]]
        gc.collect()
        # 1. create slides from coordinates
        # 2. scale to original area
        # 3. generate perineurium
        # 4. don't need to update epineurium border since this will get done later
        defrefsection = deepcopy(defsection)
        defmap = FascicleConnectivityMap(defsection)
        defmap.generate_map()
        refmap = defmap.get_simplemap()
        defdir = os.path.join(root, params['path']['slides'], 'doubledeform', 'postdeformtxts')
        deformslides = []
        for i in range(1, len(defsection) - 1):
            fascdir = os.path.join(defdir, 'p', str(i))
            nervedir = os.path.join(defdir, 'n', str(i))
            fascicles = []
            for file in [f for f in os.listdir(fascdir) if f.endswith('.txt')]:
                with open(os.path.join(fascdir, file)) as f:
                    # load in python list
                    fascpts = eval(f.read())
                # append 0 to each point
                fascpts = [[pt[0], -pt[1], 0] for pt in fascpts]
                if Trace(fascpts).area() > 1000 / 1000 / 1000:  # ignores little holes in fascicles touching
                    fascicles.append(Fascicle(Trace(fascpts)))
                else:
                    print('Drop')
            nervefiles = [f for f in os.listdir(nervedir) if f.endswith('.txt')]
            assert len(nervefiles) == 1
            with open(os.path.join(nervedir, nervefiles[0])) as f:
                # load in python list
                nervepts = eval(f.read())
            # append 0 to each point
            nervepts = [[pt[0], -pt[1], 0] for pt in nervepts]
            nerve = Nerve(Trace(nervepts))
            slide: Slide = Slide(
                fascicles,
                nerve,
                nerve_mode,
                will_reposition=(deform_mode != DeformationMode.NONE),
            )
            slide.scale(1000)
            # shift slide about (0,0)
            slide.smooth_traces(params["preprocess"]["nsmoothing"], params["preprocess"]["fsmoothing"])
            slide.move_center(np.array([0, 0]))
            deformslides.append(slide)
            assert np.all(np.isclose(slide.nerve.centroid(), [0, 0]))
        # %%
        for i, slidetwo in enumerate(zip(deformslides, defrefsection[1:-1])):
            try:
                assert len(slidetwo[0].fascicles) == len(slidetwo[1].fascicles)
            except AssertionError:
                plt.title(f'deformed: {i}')
                slidetwo[0].plot()
                plt.title(f'reference {i}')
                slidetwo[1].plot()
                test = watershed_slide(slidetwo[0], None, slidetwo[1], blind=True, plot=True, interactive=True)
                assert len(test.fascicles) == len(slidetwo[1].fascicles)
                deformslides[i] = test

        # %%
        gc.collect()
        outslides = deform_postfix_area_only(deepcopy(deformslides), defrefsection[1:-1], plot=True)
        slides[deformspan1[0] + 1 : deformspan2[1] - 1] = outslides
    slidesave = deepcopy(slides)
    # %%
    slides = deepcopy(slidesave)
    gc.collect()
    if deform_mode == DeformationMode.PHYSICS and rundeform:
        slidemods = []
        # temp
        bottom = math.floor(slidespan1[0]) - transition_n
        bottom = bottom if bottom > 0 else 0
        slidemods.append(slides[bottom : math.floor(slidespan1[0])])
        radii = [r1]
        slidemods[0].reverse()
        if cuff_config.get('zspan2') is not None:
            top = math.ceil(slidespan2[1]) + transition_n
            top = top if top < len(slides) else len(slides) - 1
            slidemods.append(slides[math.ceil(slidespan2[1]) : top])
            radii.append(r2)
        if 'adaptive' in model_config['cuff']['preset']:
            minrad = min(radii)
            if minrad + 100 < 1500:
                cuffnum = 2000
            else:
                cuffnum = 3000
            model_config["cuff"]["preset"] = f"LivaNova{cuffnum}_v2_3D.json"
            model_config["ascent_cuff"] = f"LivaNova{cuffnum}_v2.json"

        for d in slides[slidespan1[0] : slidespan1[1]]:
            d.nerve = d.reshaped_nerve(ReshapeNerveMode.CIRCLE, override_r=r1)
            d.nerve.shift([-d.nerve.centroid()[0], -d.nerve.centroid()[1], 0])
            assert np.all(np.isclose(d.nerve.centroid(), [0, 0]))
        if cuff_config.get('zspan2') is not None:
            for d in slides[slidespan2[0] : slidespan2[1]]:
                d.nerve = d.reshaped_nerve(ReshapeNerveMode.CIRCLE, override_r=r2)
                d.nerve.shift([-d.nerve.centroid()[0], -d.nerve.centroid()[1], 0])
                assert np.all(np.isclose(d.nerve.centroid(), [0, 0]))
        print('Deforming cuff ends')
        for r, ss in zip(radii, slidemods):
            slidelist = util.deforms_linear(ss, r)
        print('Deforming between the cuff centers')  # todo - adapt linear deforms to be able to do this
        if cuff_config.get('zspan2') is not None:
            rlist = []
            for i in range(len(slides[slidespan1[1] : slidespan2[0]])):
                distance = (i + 1) / len(slides[slidespan1[1] : slidespan2[0]])
                target_r = r2 * distance + r1 * (1 - distance)
                rlist.append(target_r)
            num = int((slidespan2[0] - slidespan1[1]) / 2)
            final_defrat = 1 - transition_n / (slidespan2[0] - slidespan1[1]) / 2
            slidetest = util.deforms_linear(
                slides[slidespan1[1] : slidespan1[1] + num], r1, end_deform_ratio=final_defrat, rlist=rlist[:num]
            )
            util.deforms_linear(
                list(reversed(slides[slidespan2[0] - num : slidespan2[0]])),
                r2,
                end_deform_ratio=final_defrat,
                rlist=list(reversed(rlist[-num:])),
            )
        issues = 0
        print('done deforming')
        for i, slide in enumerate(slides):
            if not slide.validate(die=False, plot_debug=True, intersection_target='inners'):
                issues += 1
                print(i)
        # %%
        # TODO: generate video from this
        if issues > 0:
            raise RuntimeError("issues found")
        # TODO: rework to use new fascicle shape. smooth and then scale to original area
    # %% Compute parameters and generate comsol file
    os.chdir(os.path.split(scriptloc)[0])
    model_config["nerve_length"] = (len(slides)) * params["input"]["um_per_slice"]
    # get medium params and save
    model_config = compute_medium_params(model_config, len(slides), params["input"]["um_per_slice"])
    with open(params['project_path'] + '/model.json', 'w') as f:
        json.dump(model_config, f, indent=2)
    if not skipsave:
        # save slides
        print("Saving slides...")
        with open(os.path.join(params['path']['slides'], 'slides.obj'), 'wb') as f:
            pickle.dump(slides, f)

        # save images now that deformation is complete
        print("Generating output images...")
        final_fmap = FascicleConnectivityMap(slides)
        final_fmap.save_images(  # TODO: make this a static method
            params['path']['slides'],
            dims,
            separate=True,
            print_ids=False,
            resize_factor=1 / params["output"]["image_um_per_px"],
            id_font=params["output"]["id_font"],
        )
        # fill holes
        for dire in ['p', 'n', 'i']:
            fill_img_holes(os.path.join(params['path']['slides'], dire))
# TODO: add fill holes to image export
# %%
if geometry:
    # temporary measure. If deforming and using imthera cuff, then open model config and edit min_radius_enclosing_circle1 to be the maximum between the two radii
    with open(os.path.join(params['project_path'], 'model.json')) as f:
        model_config = json.load(f)
    if 'ImThera' in model_config['cuff']['preset']:
        print("Temporarily editing model config to use max radius")
        model_config['min_radius_enclosing_circle1'] = max(
            model_config['min_radius_enclosing_circle1'], model_config['min_radius_enclosing_circle2']
        )
        with open(os.path.join(params['project_path'], 'model.json'), 'w') as f:
            json.dump(model_config, f, indent=2)
    print("Generating geometry in comsol...")
    run_comsol(config3d_path, "geometry")
    # load and dump matmap.json with indent=2
    with open(os.path.join(params['path']['comsol'], 'stl', 'matmap.json')) as f:
        matmap = json.load(f)
    with open(os.path.join(params['path']['comsol'], 'stl', 'matmap.json'), 'w') as f:
        json.dump(matmap, f, indent=2)
    with open(os.path.join(params['path']['comsol'], 'stl', 'dommap.json')) as f:
        dommap = json.load(f)
    # dommap is a list of objects and their domains. Remove domains which are repeated, unless they are the first domain
    removed = True
    while removed:
        removed = False
        for name, doms in dommap.items():
            if len(doms) == 1:
                continue
            for name2, doms2 in dommap.items():
                if name != name2 and len(doms2) == 1:
                    if doms2[0] in doms:
                        doms.remove(doms2[0])
                        removed = True
    with open(os.path.join(params['path']['comsol'], 'stl', 'dommap.json'), 'w') as f:
        json.dump(dommap, f, indent=2)
# %% temp check for cliffs
if slidegen and not skipsave:
    params['openbuffer'] = 10
    openbuffer = params['openbuffer']  # distance in um that peri will be shrun kwhen checking for exposed endo
    w = (dims[0][1] - dims[0][0]) / params["output"]["image_um_per_px"]
    h = (dims[1][1] - dims[1][0]) / params["output"]["image_um_per_px"]
    capvol = np.zeros((int(h), int(w), len(slides)))

    def openanalyze(a, b, openbuffer=0):
        bpers = []
        for f in b.fascicles:
            thisf = f.outer.deepcopy()
            thisf.offset(distance=-openbuffer)
            bpers.append(thisf.polygon())
        try:
            bper = unary_union(bpers)
        except ValueError:
            from shapely import validation

            bpers = [validation.make_valid(p) for p in bpers]
            bper = unary_union(bpers)
        exposed = []
        for f in a.fascicles:
            fascgood = True
            for i in f.inners:
                if not bper.contains(i.polygon()):
                    fascgood = False
            if not fascgood:
                exposed.append(f)
        return exposed

    def opencheck(a, b, openbuffer=0):
        bpers = []
        for f in b.fascicles:
            thisf = f.outer.deepcopy()
            thisf.offset(distance=-openbuffer)
            bpers.append(thisf.polygon())
        bper = unary_union(bpers)
        aendo = unary_union([i.polygon() for f in a.fascicles for i in f.inners])
        return bper.contains(aendo)

    def add_fasc_caps(capvol, issue, flist, i):
        imbuffer = 0

        def prep_points(points):
            # adjusts plot points to dimensions and formats for PIL
            points = (points - dim_min + imbuffer)[:, 0:2].astype(int)
            points = tuple(zip(points[:, 0], points[:, 1]))
            return points

        def factor_resize(im, factor):
            (width, height) = (int(im.width * factor), int(im.height * factor))
            im_resized = im.resize((width, height))
            return im_resized

        from PIL import Image, ImageDraw

        resize_factor = 1 / params["output"]["image_um_per_px"]
        um_per_slice = params['input']["um_per_slice"]
        dim_min = [min(x) for x in dims]
        dim = [max(x) for x in dims]
        imdim = [dim[0] + abs(dim_min[0]) + imbuffer * 2, dim[1] + abs(dim_min[1]) + imbuffer * 2]
        for f in flist:
            added_thickness = 0
            added_slices = 0
            while added_thickness < f.outer.thickness:
                imgp = Image.new('L', imdim)
                draw = ImageDraw.Draw(imgp)
                pixel_brightness = 256 * (f.outer.thickness - added_thickness) / um_per_slice
                value = 256 if pixel_brightness > 256 else pixel_brightness
                value = 128 if value < 128 else value
                added_thickness += um_per_slice * value / 256
                this_slice = i + issue + issue * added_slices
                added_slices += 1
                draw.polygon(prep_points(f.outer.points[:, 0:2]), fill=int(value))
                imgp = imgp.transpose(Image.Transpose.FLIP_TOP_BOTTOM)
                imgp = factor_resize(imgp, resize_factor)
                imadd = np.array(imgp)
                capvol[:, :, this_slice] += imadd

    badtransition = []
    for i in range(len(slides)):
        if i > 0:
            exposed_fascicles = openanalyze(slides[i], slides[i - 1], openbuffer)
            if len(exposed_fascicles) > 0:
                badtransition.append(-i)
                issue = -1
                add_fasc_caps(capvol, issue, exposed_fascicles, i)
        if i < len(slides) - 1:
            exposed_fascicles = openanalyze(slides[i], slides[i + 1], openbuffer)
            if len(exposed_fascicles) > 0:
                badtransition.append(i)
                issue = 1
                add_fasc_caps(capvol, issue, exposed_fascicles, i)
    pd.DataFrame(badtransition).to_csv(params['path']['mesh'] + '/badtransition.csv', index=False)
    outpath = params['path']['slides'] + '/pcap'
    ensure_dir(outpath)
    for i in range(len(slides)):
        im = capvol[:, :, i].astype(int)
        cv2.imwrite(outpath + f'/{i}.tif', im)

# %% Set up config for SIP premesh
print("Running mesh in ScanIP...")

with open(root + '/config/run.json') as f:
    run_config = json.load(f)
if mesh:
    meshstart = time.time()
    sipconfig = {
        'mesh': params['mesh'],
        'outpath': params['path']['mesh'],
        'n_imgs': params['path']['slides'] + '/n',
        'i_imgs': params['path']['slides'] + '/i',
        'p_imgs': params['path']['slides'] + '/p',
        'pcap_imgs': params['path']['slides'] + '/pcap',
        'um_per_px': params["output"]["image_um_per_px"],
        'um_per_slice': params['input']['um_per_slice'],
        'run_type': run_config['submission_context'] if quit_premesh is not True else "local",
        'stl_path': params['path']['comsol'] + '/stl',
        'sep_nerve': ascent_configs.search(Config.SAMPLE, "boundary_separation", "nerve"),
    }
    sip_config_path = root + '/config/system/sipmeshconfig.json'
    scriptname = 'sip_mesher_ds'
    run_scanip(sip_config_path, sipconfig, scriptname, progress=args.mesher_progress)
    meshend = time.time()
    mesh_hrs = (meshend - meshstart) / 60 / 60
    np.savetxt(params['path']['mesh'] + '/mesh_hours.dat', [mesh_hrs])

if quit_premesh:
    sys.exit('Quitting premesh')

# %%
if model and run_config['submission_context'] == "cluster":
    print('Generating Comsol model.')
    # prior to running java model code, goes through nastran mesh and gets material/contact ids
    meshfile = params['path']['mesh'] + '/mesh.nas'
    contact_ids = {}
    material_ids = {}
    with open(meshfile) as m:
        meshdata = list(m.readlines())
    for i, line in enumerate(meshdata):
        if '$Material Index:' in line:
            index = line.split()[-1]
            matline = meshdata[i - 2]
            mat = matline.split()[-1]
            if material_ids.get(mat) is None:
                material_ids[mat] = [index]
            else:
                material_ids[mat].append(index)
        if 'Shell definition for contact surface Contact' in line:
            contact = line.split(':')[-1][1:-1]
            index = meshdata[i + 1].split(',')[1]
            contact_ids[contact] = index
    with open(params['path']['comsol'] + '/stl/meshmap.json', 'w') as out:
        json.dump({'materials': material_ids, 'contacts': contact_ids}, out, indent=2)
    solvestart = time.time()
    run_comsol(config3d_path, "model")
    solveend = time.time()
    solve_hours = (solveend - solvestart) / 60 / 60
    np.savetxt(params['path']['comsol'] + '/solve_hours.dat', [solve_hours])
# %%Set up ascent directory and prep for fibersets
print('Generating files for ascent')
if os.path.exists(params['path']['ascent']):
    shutil.rmtree(params['path']['ascent'])
os.chdir(root)
configdir = os.path.split(params['config']['path'])[1]
zlocs = params['output']['slices']
slices = zlocs.keys()

ensure_dir(params['path']['fibers'] + '/3D_fiberset')
ensure_dir(params['path']['ascent'] + '/inputs')
ensure_dir(params['path']['ascent'] + '/samples')
ensure_dir(params['path']['ascent'] + '/samples/3D/models/0/sims/3')
ensure_dir(params['path']['ascent'] + '/samples/3D/models/0/sims/3/ss_bases')
ensure_dir(params['path']['ascent'] + '/samples/3D/models/0/sims/3/ss_coords')
ensure_dir(params['path']['ascent'] + '/samples/3D/models/0/sims/3/ss_lengths')

shutil.copyfile(
    os.path.join(params['project_path'], 'model.json'),
    f"{params['path']['ascent']}/samples/3D/models/0/model.json",
)
shutil.copyfile(
    f"{configdir}/{samplefile}",
    f"{params['path']['ascent']}/samples/3D/sample.json",
)
zs = []
# get pcs z locations
for sli in slices:
    this_pseudo = pseudonym + '-' + sli
    ensure_dir(params['path']['ascent'] + '/inputs/' + this_pseudo)
    ensure_dir(params['path']['ascent'] + '/samples/' + this_pseudo)
    ensure_dir(params['path']['ascent'] + '/samples/' + this_pseudo + '/models/0/sims/3')
    # update model con
    with open(os.path.join(params['project_path'], 'model.json')) as m:
        model_config = json.load(m)
    # read in min_radius set during slidegen step and set as r_nerve_override
    # TODO: if deformation set this to mean area under cuff used for deformation?
    # TODO make this more robust, but for now passing
    if False:
        if params['input']['z-'] in sli:
            if deform_mode == DeformationMode.PHYSICS:
                model_config['3d_nerve_override'] = r1
            else:
                model_config['3d_nerve_override'] = model_config['min_radius_enclosing_circle1']
        else:
            if deform_mode == DeformationMode.PHYSICS:
                model_config['3d_nerve_override'] = r2
            else:
                model_config['3d_nerve_override'] = model_config['min_radius_enclosing_circle2']
    else:
        print(
            'NOTEWARNINGNOTE: skipping output of 3D nerve override, must manually put this value in model when running paired ASCENT models'
        )
    with open(f"{params['path']['ascent']}/samples/{this_pseudo}/models/0/model.json", 'w') as m:
        model_config["cuff"]["preset"] = model_config["ascent_cuff"]
        json.dump(model_config, m, indent=2)
    # update sample con
    with open(f"{configdir}/{samplefile}") as m:
        sample_config = json.load(m)
    sample_config["boundary_separation"]["nerve"] = 0
    sample_config["boundary_separation"]["fascicles"] = 0
    sample_config["sample"] = this_pseudo
    sample_config["scale"]["shrinkage"] = 0
    with open(f"{params['path']['ascent']}/samples/{this_pseudo}/sample.json", 'w') as m:
        json.dump(sample_config, m, indent=2)

    for file in os.listdir(params['path']['comsol'] + '/pcs'):
        if file.startswith('pcs'):
            zs.append(np.loadtxt(params['path']['comsol'] + '/pcs/' + file)[-1])
with open(params['path']['fibers'] + '/zlocs.json', 'w') as f:
    json.dump(zlocs, f, indent=2)
zslices = {}
# get slides
for sli, zval in zlocs.items():
    this_pseudo = pseudonym + '-' + sli
    index = int(round(zval / params["input"]["um_per_slice"], 0))
    zslices[sli] = index
    for x in ['i', 'n', 'p']:
        shutil.copyfile(
            params['path']['slides'] + f'/{x}/{index}.tif',
            params['path']['ascent'] + '/inputs/' + this_pseudo + f"/{x if x != 'p' else 'o'}.tif",
        )
with open(params['path']['fibers'] + '/zslices.json', 'w') as f:
    json.dump(zslices, f, indent=2)
shutil.copyfile(params['path']['fibers'] + '/zslices.json', params['path']['ascent'] + '/zslices.json')
shutil.copyfile(params['path']['fibers'] + '/zlocs.json', params['path']['ascent'] + '/zlocs.json')
# %% generate fiberset and extract all potentials


def get_pots(mphfile, fiberspath, outpath):
    # run code for potential extraction
    import subprocess

    comsol_path = env['ASCENT_COMSOL_PATH']
    jdk_path = env['ASCENT_JDK_PATH']
    project_path = root + '/src/comsol_java_api'
    core_name = 'PotentialsExtractor'

    if sys.platform.startswith('darwin'):  # macOS
        subprocess.Popen([f'{comsol_path}/bin/comsol', 'server'], close_fds=True)
        time.sleep(5)
        os.chdir(project_path + '/' + 'src')
        os.system(
            f'{jdk_path}/javac -classpath ../bin/json-20190722.jar:{comsol_path}/plugins/* model/*.java -d ../bin'
        )
        # https://stackoverflow.com/questions/219585/including-all-the-jars-in-a-directory-within-the-java-classpath
        os.system(
            f'{comsol_path}/java/maci64/jre/Contents/Home/bin/java '
            f'-cp .:$(echo {comsol_path}/plugins/*.jar | '
            f'tr \' \' \':\'):../bin/json-20190722.jar:../bin model.{core_name} "{mphfile}" "{fiberspath}" "{outpath}"'
        )
        os.chdir('..')

    elif sys.platform.startswith('linux'):  # linux
        subprocess.Popen([f'{comsol_path}/bin/comsol', 'server'], close_fds=True)
        time.sleep(5)
        os.chdir(project_path + '/' + 'src')
        print(os.getcwd())
        os.system(
            f'{jdk_path}/javac -classpath ../bin/json-20190722.jar:{comsol_path}/plugins/* model/*.java -d ../bin'
        )
        # https://stackoverflow.com/questions/219585/including-all-the-jars-in-a-directory-within-the-java-classpath
        os.system(
            f'{comsol_path}/java/glnxa64/jre/bin/java '
            f'-cp .:$(echo {comsol_path}/plugins/*.jar | '
            f'tr \' \' \':\'):../bin/json-20190722.jar:../bin model.{core_name} "{mphfile}" "{fiberspath}" "{outpath}"'
        )
        os.chdir('..')

    else:  # assume to be 'win64'
        subprocess.Popen([f'{comsol_path}\\bin\\win64\\comsolmphserver.exe'], close_fds=True)
        time.sleep(5)
        os.chdir(project_path + '/' + 'src')
        os.system(
            f'""{jdk_path}\\javac" '
            f'-cp "..\\bin\\json-20190722.jar";"{comsol_path}\\plugins\\*" '
            f'model\\*.java -d ..\\bin"'
        )
        os.system(
            f'""{comsol_path}\\java\\win64\\jre\\bin\\java" '
            f'-cp "{comsol_path}\\plugins\\*";"..\\bin\\json-20190722.jar";"..\\bin" '
            f'model.{core_name} "{mphfile}" "{fiberspath}" "{outpath}""'
        )
        os.chdir('..')


if fibergen:
    print('Generating fiberset')
    run_comsol(config3d_path, "fibers")
if extract:
    print('Fiberset generation complete')
    fibertxt = params['path']['fibers'] + '/fibers.txt'
    data = pd.read_csv(fibertxt, sep=",", header=0, skiprows=7)
    data.rename(columns={'% x': 'x', 'Streamline': 'line'}, inplace=True)
    fibers = {}
    # separate fiber data into dictionary items
    for line in pd.unique(data.line):
        fiberdata = data.iloc[np.where(data.line == line)[0], :]
        if fiberdata.iloc[0]['z'] > 100:
            continue
        fibers[line] = fiberdata
    print(f'Generating {len(fibers.keys())} supersampled fibers to extract from COMSOL: ')
    # resample each fiber according to sim settings
    for i, line in enumerate(fibers.keys()):
        print(i, end=' ')
        fiberdata = fibers[line]
        fiberpoints = np.array(fiberdata.loc[:, ['x', 'y', 'z']])
        fiberflip = np.flip(fiberpoints, axis=0)
        nd = nd_line(fiberflip)
        le = nd.length
        ascent_configs.configs['sims']['fibers']['z_parameters']['max'] = le
        ascent_configs.configs['models']['medium']['proximal']['length'] = le
        ascent_configs.configs['models']['medium']['proximal']['length'] = le
        np.savetxt(params['path']['ascent'] + f'/samples/3D/models/0/sims/3/ss_lengths/{i}.dat', [le])
        fiberset = FiberSet(None)
        fiberset.configs.update(ascent_configs.configs)
        fpoints = fiberset._generate_z([(0, 0)], super_sample=True, override_length=le)
        fpoints = [d['fiber'] for d in fpoints]
        fpoints = np.vstack(fpoints)[:, -1]
        coords = np.zeros([len(fpoints), 3])
        coords[:, 2] = fpoints
        np.savetxt(
            params['path']['ascent'] + f'/samples/3D/models/0/sims/3/ss_coords/{i}.dat',
            coords,
            delimiter=' ',
            fmt='%.10f',
            header=str(len(fpoints)),
            comments='',
        )
        sample_points = [nd.interp(p) for p in fpoints]
        # END TEMP
        np.savetxt(
            params['path']['fibers'] + f'/3D_fiberset/{i}.dat',
            sample_points,
            delimiter=' ',
            fmt='%.10f',
            header=str(len(sample_points)),
            comments='',
        )
    else:
        print('')
    # format fiber xy points for ascent
    for file in slices:
        this_pseudo = pseudonym + '-' + file
        fibertxt = params['path']['fibers'] + f'/{file}.txt'
        data = pd.read_csv(fibertxt, sep=",", header=0, skiprows=7)
        outfib = []
        for i, row in data.iterrows():
            if i in fibers:
                outfib.append(np.array([row['% x start'], row['y start']]))
        ensure_dir(params['path']['ascent'] + f'/samples/{this_pseudo}/explicit_fibersets')
        np.savetxt(
            params['path']['ascent'] + f'/samples/{this_pseudo}/explicit_fibersets/0.txt',
            outfib,
            delimiter=' ',
            fmt='%.10f',
            header='%% Fiber XY Coordinates: X (space) Y',
            comments='',
        )
    # extract potentials
    pcspath = params['path']['comsol'] + '/pcs'
    numpcs = len([f.split('.')[0] for f in os.listdir(pcspath) if f.endswith('.dat')])
    for i in range(numpcs):
        mph = params['path']['comsol'] + f'/{i}.mph'
        fibers = params['path']['fibers'] + '/3D_fiberset'
        out = params['path']['comsol'] + f'/potentials/{i}'
        ensure_dir(out)
        if not skip_getpots:
            get_pots(mph, fibers, out)
        if os.path.exists(params['path']['ascent'] + f'/samples/3D/models/0/sims/3/ss_bases/{i}'):
            shutil.rmtree(params['path']['ascent'] + f'/samples/3D/models/0/sims/3/ss_bases/{i}')
        shutil.copytree(out, params['path']['ascent'] + f'/samples/3D/models/0/sims/3/ss_bases/{i}')
    shutil.copytree(fibers, params['path']['ascent'] + '/samples/3D/models/0/sims/3/3D_fiberset')
    ensure_dir(env['THREED_DATA_EXPORT_PATH'])
    # copy ascent output from ascent path to out
    if os.path.exists(env['THREED_DATA_EXPORT_PATH'] + f'/{pseudonym}'):
        shutil.rmtree(env['THREED_DATA_EXPORT_PATH'] + f'/{pseudonym}')
    shutil.copytree(params['path']['ascent'], env['THREED_DATA_EXPORT_PATH'] + f'/{pseudonym}')
# %%
print('fin')
