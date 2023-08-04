# # -*- coding: utf-8 -*-
# """
# Created on Tue Oct 19 14:09:22 2021

# @author: dpm42
# """
# import os

# os.chdir('d:/ASCENt/fresh')

import sys

c = input(
    '''
          Which config would you like to edit?
          r=Run
          s=Sim
          e=Sample
          m=Model
          '''
)

import json


def param_picker(file):
    with open(file, 'r') as f:
        config = json.load(f)
    choice = [input(f'Pick a parameter to change (or a subparam) from the list:\n{config.keys()}')]
    while type(config[choice[-1]]) == dict:
        config = config[choice[-1]]
        choice.append(input(f'Pick a parameter to change (or a subparam) from the list:\n{config.keys()}'))
    print(f'Current param value for {choice}:{config[choice[-1]]}')
    newval = input('New value for param:')
    return choice, newval


def configval(path, param, val):
    with open(path, 'r') as f:
        config = json.load(f)
    cf = config
    for p in param[:-1]:
        cf = cf.get(p)
    if cf is None or cf.get(param[-1]) is None:
        sys.exit('Invalid param')
    oldval = cf.get(param[-1])
    cf.update({param[-1]: val})
    with open(path, 'w') as f:
        json.dump(config, f, indent=2)
    return oldval


import os
import shutil

os.chdir('..')

if c == 'm':
    samples =[
    370,
    3701,
    3719,
    3721,
    650]
    models = [0]

    change = ['inner_interp_tol']

    newval = 0.015

    clear = False

    for sample in samples:
        for model in models:
            if clear:
                try:
                    os.remove(f'samples/{sample}/models/{model}/debug_geom.mph')
                    shutil.rmtree(f'samples/{sample}/models/{model}/sims')
                except:
                    print(f'Could not delete for {sample}, {model}.')
            try:
                oldval = configval(f'samples/{sample}/models/{model}/model.json', change, newval)
                print(f'Updated sample {sample} model {model} parameter {change} from {oldval} to {newval}.')
            except Exception as e:
                print(sample, e)

#%%
elif c == 's':
    sims = [1, 10, 11, 12, 13, 14]

    change = ["fibers", "xy_parameters", "r_shift"]

    newval = 2390.963

    for sim in sims:
        try:
            oldval = configval(f'config/user/sims/{sim}.json', change, newval)
            print(f'Updated sim {sim} parameter {change} from {oldval} to {newval}.')
        except:
            pass
#%%
elif c == 'r':
    path = 'config/user/runs'

    # runlist = [os.path.splitext(x)[0] for x in os.listdir(path) if x.endswith('.json')]

    # run = input(f'Pick from the available runs:\n{runlist}')

    runlist = [
        '250',
        '2501',
        '2509',
        '251',
        '252',
        '2521',
        '2529',
        '253',
        '2531',
        '270',
        '271',
        '272',
        '273',
        '370',
        '3701',
        '3709',
        '371',
        '372',
        '3721',
        '3729',
        '373',
        '3731',
        '570',
        '5701',
        '5709',
        '571',
        '572',
        '5721',
        '5729',
        '573',
        '5731',
        '650',
        '651',
        '652',
        '653',
        '670',
        '6701',
        '6709',
        '671',
        '672',
        '6721',
        '6729',
        '673',
        '6731',
        '2519',
        '2511',
        '3711',
        '3719',
        '6711',
        '6719',
        '5711',
        '5719',
    ]

    change = ['sims']

    # newval = [3]

    print(str([int(x) for x in runlist]))

    print(str([int(x) for x in runlist]).replace(',', ''))

    # newval = [3, 7, 10,11,20,30]
    newval = [3, 10]


    for run in runlist:
        # change, newval = param_picker(path + f'/{run}.json')

        oldval = configval(f'config/user/runs/{run}.json', change, newval)

        print(f'Updated run {run} parameter {change} from {oldval} to {newval}.')
#%%
elif c == 'e':
    samples = [250,
     2501,
     2509,
     251,
     2511,
     2515,
     2519,
     252,
     2520,
     2521,
     2524,
     2526,
     2529,
     253,
     2530,
     2531,
     2534,
     2535,
     2536,
     270,
     271,
     272,
     273,
     370,
     3701,
     3709,
     371,
     3711,
     3719,
     372,
     3721,
     3729,
     373,
     3731,
     570,
     5701,
     5709,
     571,
     5711,
     5715,
     5719,
     572,
     5721,
     5729,
     573,
     5731,
     5735,
     650,
     651,
     652,
     653,
     670,
     6701,
     6709,
     671,
     6711,
     6719,
     672,
     6721,
     6729,
     673,
     6731]
    change = ["boundary_separation"]

    newval =  {'fascicles':10,'nerve':10}

    clear = False

    for sample in samples:
        try:
            oldval = configval(f'samples/{sample}/sample.json', change, newval)
            print(f'Updated sample {sample} parameter {change} from {oldval} to {newval}.')
        except:
            print('oopsie')

else:
    print('Not implemented')
