#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

# builtins
import os
import sys
import time

# ascent
from src.runner import Runner
from src.utils.enums import Config, Env, SetupMode

from .env_setup import run as env_setup


def run(args):
    # test
    if not (sys.version_info.major == 3 and sys.version_info.minor >= 7):
        print(f'You are running Python {sys.version_info.major}.{sys.version_info.minor}, but 3.7 or later required')
        sys.exit()

    # create bin/ directory for storing compiled Java files if it does not yet exist
    if not (os.path.exists('bin')):
        os.mkdir('bin')

    for argument in args.run_indices:
        # START timer
        start = time.time()

        try:
            int(argument)
        except ValueError:
            print(f'Invalid type for argument: {argument}\nAll arguments must be positive integers.')
            sys.exit()

        if int(argument) < 0:
            print(f'Invalid sign for argument: {argument}\nAll arguments must be positive integers.')
            sys.exit()

        print(f'\n\n########## SUBMITTING RUN {argument} ##########\n\n')

        run_path = os.path.join('config', 'user', 'runs', f'{argument}.json')
        if not os.path.exists(run_path):
            print(f'Invalid run configuration path: {run_path}')
            sys.exit()

        env_path = os.path.join('config', 'system', 'env.json')
        if not os.path.exists(env_path):
            print(f'Missing env configuration file: {env_path}')
            env_setup(env_path)

        # initialize Runner (loads in parameters)
        runner = Runner(int(argument))
        runner.add(SetupMode.NEW, Config.RUN, run_path)
        runner.add(SetupMode.NEW, Config.ENV, env_path)
        runner.add(SetupMode.OLD, Config.CLI_ARGS, vars(args))

        # populate environment variables
        runner.populate_env_vars()

        # ready, set, GO!
        runner.submit(argument)

        # END timer
        end = time.time()
        elapsed = end - start

        print(f"\nruntime: {time.strftime('%H:%M:%S', time.gmtime(elapsed))} (hh:mm:ss)")

    # cleanup for console viewing/inspecting
    del start, end
