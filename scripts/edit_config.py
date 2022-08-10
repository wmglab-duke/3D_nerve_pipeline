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

os.chdir('d:/ASCENt/fresh')

if c == 'm':
    samples = [2110]
    models = [11, 12, 21, 22, 23, 31, 32, 33]

    change = ['medium', 'distal', 'radius']

    newval = 20040.963

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
            except:
                pass

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

    runlist = [os.path.splitext(x)[0] for x in os.listdir(path) if x.endswith('.json')]

    run = input(f'Pick from the available runs:\n{runlist}')

    change, newval = param_picker(path + f'/{run}.json')

    oldval = configval(f'config/user/runs/{run}.json', change, newval)

    print(f'Updated run {run} parameter {change} from {oldval} to {newval}.')

else:
    print('Not implemented')
