#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
1"""

import argparse
import json
import multiprocessing
import os
import pickle
import re
import shutil
import subprocess
import sys
import time
import warnings

import numpy as np
import pandas as pd


# %%Set up parser and top level args
class listAction(argparse.Action):
    def __call__(self, parser, values, option_string=None):
        run_path = 'runs'
        jsons = [file for file in os.listdir(run_path) if file.endswith('.json')]
        data = []
        for j in jsons:
            with open(run_path + '/' + j) as f:
                try:
                    rundata = json.load(f)
                except Exception as e:
                    print('WARNING: Could not load {}'.format(j))
                    print(e)
                    continue
                data.append(
                    {
                        'RUN': os.path.splitext(j)[0],
                        'PSEUDONYM': rundata.get('pseudonym'),
                        'SAMPLE': rundata['sample'],
                        'MODELS': rundata['models'],
                        'SIMS': rundata['sims'],
                    }
                )
        df = pd.DataFrame(data)
        df.RUN = df.RUN.astype(int)
        df = df.sort_values('RUN')
        print('Run indices available (defined by user .json files in {}):\n'.format(run_path))
        print(df.to_string(index=False))
        sys.exit()


parser = argparse.ArgumentParser(
    description='ASCENT: Automated Simulations to Characterize Electrical Nerve Thresholds'
)
parser.add_argument(
    'run_indices',
    type=int,
    nargs='*',
    help='Space separated indices to submit NEURON sims for',
)
parser.add_argument('-p', '--partition', help='If submitting on a cluster, overrides slurm_params.json')
parser.add_argument(
    '-n',
    '--num-cpu',
    type=int,
    help='For local submission: set number of CPUs to use, overrides run.json',
)
parser.add_argument(
    '-m',
    '--job-mem',
    type=int,
    help='For cluster submission: set amount of RAM per job (in MB), overrides slurm_params.json',
)
parser.add_argument(
    '-j',
    '--num-jobs',
    type=int,
    help='For cluster submission: set number of jobs per array, overrides slurm_params.json',
)
parser.add_argument(
    '-l',
    '--list-runs',
    action=listAction,
    nargs=0,
    help='List info for available runs.z If supplying this argument, do not pass any run indices',
)
parser.add_argument(
    '-A',
    '--all-runs',
    action='store_true',
    help='Submit all runs in the present export folder. If supplying this argument, do not pass any run indices',
)
parser.add_argument(
    '-s',
    '--skip-summary',
    action='store_true',
    help='Begin submitting fibers without asking for confirmation',
)
parser.add_argument(
    '-S',
    '--slurm-params',
    type=str,
    help='For cluster submission: string for additional slurm parameters (enclose in quotes)',
)
parser.add_argument(
    '-c',
    '--force-recompile',
    action='store_true',
    help='Force submit.py to recompile NEURON files',
)
submit_context_group = parser.add_mutually_exclusive_group()
submit_context_group.add_argument(
    '-L',
    '--local-submit',
    action='store_true',
    help='Set submission context to local, overrides run.json',
)
submit_context_group.add_argument(
    '-C',
    '--cluster-submit',
    action='store_true',
    help='Set submission context to cluster, overrides run.json',
)

parser.add_argument('-v', '--verbose', action='store_true', help='Print detailed submission info')

OS = 'UNIX-LIKE' if any([s in sys.platform for s in ['darwin', 'linux']]) else 'WINDOWS'


# %% Set up utility functions


class WarnOnlyOnce:
    warnings = set()

    @classmethod
    def warn(cls, message):
        # storing int == less memory then storing raw message
        h = hash(message)
        if h not in cls.warnings:
            # do your warning
            print(f"Warning: {message}")
            cls.warnings.add(h)


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█'):
    """Print or update a progress bar in the terminal.
    Call in a loop to create a terminal progress bar.

    :param iteration: The current iteration (current/total)
    :param total: The total number of iterations
    :param prefix: The prefix string to place before the progress bar
    :param suffix: The suffix string to place after the progress bar
    :param decimals: The number of decimals to show on the percentage progress
    :param length: The length of the progress bar
    :param fill: The character to fill the progress bar with
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print('\r{} |{}| {}% {}'.format(prefix, bar, percent, suffix), end='')
    # Print New Line on Complete
    if iteration == total:
        print()


def load(config_path: str):
    """Loads in json data and returns to user, assuming it has already been validated.
    :param config_path: the string path to load up
    :return: json data (usually dict or list)
    """
    with open(config_path, "r") as handle:
        return json.load(handle)


def ensure_dir(directory):
    """Ensure that a directory exists. If it does not, create it.
    :param directory: the string path to the directory
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def auto_compile(override: bool = False):
    if (not os.path.exists(os.path.join('x86_64')) and OS == 'UNIX-LIKE') or \
            (not os.path.exists(os.path.join('MOD_Files', 'nrnmech.dll')) and OS == 'WINDOWS') or override:
        print('compile')
        os.chdir(os.path.join('MOD_Files'))
        exit_data = subprocess.run(['nrnivmodl'], shell=True, capture_output=True, text=True)
        if exit_data.returncode != 0:
            print(exit_data.stderr)
            sys.exit("Error in compiling of NEURON files. Exiting...")
        os.chdir('..')
        shutil.copytree('MOD_Files/x86_64', 'x86_64')
        compiled = True
    else:
        print('skipped compile')
        compiled = False
    return compiled


def make_task(my_os: str, sub_con: str, start_p: str, sim_p: str, fiber_path: str, inner_ind: int, fiber_ind: int,
              potentials_path: str, waveform_path: str):
    with open(start_p, 'w+') as handle:
        if my_os == 'UNIX-LIKE':
            lines = [
                '#!/bin/bash\n',
                'cd ../../\n',
                'python NEURON_Files/run_controls.py \"{}\" \"{}\" \"{}\" \"{}\" \"{}\" \"{}\"\n' \
                    .format(fiber_path,
                            inner_ind,
                            fiber_ind,
                            potentials_path,
                            waveform_path,
                            sim_p),
            ]

            if sub_con == 'cluster':
                lines.remove('cd ../../\n')

        # else:  # OS is 'WINDOWS'
        #     sim_path_win = os.path.join(*sim_p.split(os.pathsep)).replace('\\', '\\\\')
        #     lines = [
        #         'nrniv -nobanner '
        #         '-dll \"{}/MOD_Files/nrnmech.dll\" '
        #         '-c \"strdef sim_path\" '
        #         '-c \"sim_path=\\\"{}\"\" '
        #         '-c \"inner_ind={}\" '
        #         '-c \"fiber_ind={}\" '
        #         '-c \"stimamp_top={}\" '
        #         '-c \"stimamp_bottom={}\" '
        #         '-c \"fiberD={}\" '
        #         '-c \"deltaz={:.4f}\" '
        #         '-c \"axonnodes={}\" '
        #         '-c \"saveflag_end_ap_times=0\" '  # for backwards compatible, overwritten in launch.hoc if 1
        #         '-c \"saveflag_runtime=0\" '  # for backwards compatible, overwritten in launch.hoc if 1
        #         '-c \"load_file(\\\"launch.hoc\\\")\" blank.hoc\n'.format(os.getcwd(),
        #                                                                   sim_path_win,
        #                                                                   inner,
        #                                                                   fiber,
        #                                                                   top,
        #                                                                   bottom,
        #                                                                   diam,
        #                                                                   deltaz,
        #                                                                   axonnodes)
        #     ]

        handle.writelines(lines)
        handle.close()


def local_submit(fiber_data: dict):
    """Submit a fiber simulation to the local machine.
    :param fiber_data: the dictionary of fiber data for suvbmission
    """
    a = fiber_data["job_number"]
    out_path = os.path.join('logs', 'out', f'{a}.log')
    err_path = os.path.join('logs', 'err', f'{a}.log')
    start = os.path.join('start_scripts', f'start_{a}')
    with open(out_path, "w+") as fo, open(err_path, "w+") as fe:
        subprocess.run(['bash', start + '.sh'] if OS == 'UNIX-LIKE' else [start + '.bat'], stdout=fo, stderr=fe)

    # print fiber completion
    if fiber_data['verbose']:
        print(f'Completed NEURON simulation for inner {fiber_data["inner"]} fiber {fiber_data["fiber"]}.')


def submit_fibers(submission_context, submission_data):
    """Submit fiber simulations, either locally or to a cluster.
    :param submission_context: the string name of the submission_context
    :param submission_data: the dictionary of data for fiber submission
    """
    # configuration is not empty
    runcount = 0
    sim_dir = os.path.join('n_sims')

    for sim_name, runfibers in submission_data.items():
        if args.verbose:
            print('\n\n################ {} ################\n\n'.format(sim_name))

        sim_path = os.path.join(sim_dir, sim_name)
        start_dir = os.path.join(sim_path, 'start_scripts')
        start_path_base = os.path.join(start_dir, 'start_')

        if submission_context == 'cluster':

            cluster_submit(runfibers, sim_name, sim_path, start_path_base)
            runcount += 1
            if not args.verbose:
                print_progress_bar(
                    runcount, len(submission_data), length=40, prefix='Current n_sim {}:'.format(sim_name)
                )
        else:

            if args.num_cpu is not None:
                cpus = args.num_cpu

                if cpus > multiprocessing.cpu_count() - 1:
                    raise ValueError('num_cpu argument is more than cpu_count-1 CPUs')

                print(f"Submitting locally to {cpus} CPUs")

            else:
                cpus = multiprocessing.cpu_count() - 1
                warnings.warn(
                    f"You did not define number of cores to use (-n), so proceeding with cpu_core_count-1={cpus}"
                )
            os.chdir(sim_path)
            with multiprocessing.Pool(cpus) as p:
                for x in runfibers:
                    x['verbose'] = args.verbose
                if not args.verbose:
                    print_progress_bar(0, len(runfibers), length=40, prefix='n_sim {}:'.format(sim_name))
                # open pool instance, set up progress bar, and iterate over each job
                for i, _ in enumerate(p.imap_unordered(local_submit, runfibers, 1)):
                    if not args.verbose:
                        print_progress_bar(
                            i + 1,
                            len(runfibers),
                            length=40,
                            prefix='n_sim {}:'.format(sim_name),
                        )
            os.chdir("../..")


def cluster_submit(runfibers, sim_name, sim_path, start_path_base):
    """Submit fiber simulations on a slurm-based high performance computer cluster.
    :param runfibers: the list of fiber data for submission
    :param sim_name: the string name of the n_sim
    :param sim_path: the string path to the simulation
    :param start_path_base: the string prefix for all start scripts
    """
    slurm_params = load(os.path.join('config', 'system', 'slurm_params.json'))
    out_dir = os.path.abspath(os.path.join(sim_path, 'logs', 'out', ''))
    err_dir = os.path.abspath(os.path.join(sim_path, 'logs', 'err', ''))
    # assign params for array submission
    partition = slurm_params['partition'] if args.partition is None else args.partition
    njobs = slurm_params['jobs_per_array'] if args.num_jobs is None else args.num_jobs
    mem = slurm_params['memory_per_fiber'] if args.job_mem is None else args.job_mem
    array_fibertasks = [runfibers[x : x + njobs] for x in range(0, len(runfibers), njobs)]
    for tasklist in array_fibertasks:

        array_indices = [task['job_number'] for task in tasklist]

        # print fiber submission
        if args.verbose:
            for task in tasklist:
                print(f"RUNNING inner ({task['inner']}) fiber ({task['fiber']})")
                time.sleep(1)

        # submit batch job for fiber

        command = [
            'sbatch',
            *([args.slurm_params] if args.slurm_params else []),
            '--job-name={}'.format(sim_name),
            '--output={}%a.log'.format(out_dir),
            '--error={}%a.log'.format(err_dir),
            '--array={}'.format(','.join([str(x) for x in array_indices])),
            '--mem={}'.format(mem),
            '--partition={}'.format(partition),
            '--cpus-per-task=1',
            'array_launch.slurm',
            start_path_base,
        ]

        if not args.verbose:
            exit_data = subprocess.run(command, capture_output=True, text=True)
        else:
            exit_data = subprocess.run(command, capture_output=True, text=True)
            print(exit_data.stdout)
        if exit_data.returncode != 0:
            print(exit_data.stderr)
            sys.exit('Non-zero exit code during job array submission. Exiting.')

        # allow job to start before removing slurm file
        time.sleep(1.0)


def make_fiber_tasks(submission_list, submission_context):
    """Create all shell scripts for fiber submission tasks.
    :param submission_list: the list of fibers to be submitted
    """
    # assign appropriate configuration data
    sim_dir = os.path.join('n_sims')
    for sim_name, runfibers in submission_list.items():

        sim_path = os.path.join(sim_dir, sim_name)
        fibers_path = os.path.abspath(os.path.join(sim_path, 'data', 'inputs'))
        output_path = os.path.abspath(os.path.join(sim_path, 'data', 'outputs'))
        start_dir = os.path.join(sim_path, 'start_scripts')
        start_path_base = os.path.join(start_dir, 'start_')

        # ensure log directories exist
        out_dir = os.path.abspath(os.path.join(sim_path, 'logs', 'out', ''))
        err_dir = os.path.abspath(os.path.join(sim_path, 'logs', 'err', ''))
        for cur_dir in [
            fibers_path,
            output_path,
            out_dir,
            err_dir,
            start_dir,
        ]:
            ensure_dir(cur_dir)

        # load generic instance of Fiber class for given n_sim
        fiber_path = os.path.join(sim_path, 'fiber.obj')
        for fiber_data in runfibers:
            inner_ind, fiber_ind = fiber_data['inner'], fiber_data['fiber']

            start_path = '{}{}{}'.format(
                start_path_base,
                fiber_data["job_number"],
                '.sh' if OS == 'UNIX-LIKE' else '.bat',
            )

            potentials_path = os.path.join(sim_path, 'data', 'inputs',
                                           'inner{}_fiber{}.dat'.format(inner_ind, fiber_ind))
            waveform_path = os.path.join(sim_path, 'data', 'inputs', 'waveform.dat')
            make_task(OS,
                      submission_context,
                      start_path,
                      sim_path,
                      fiber_path,
                      inner_ind,
                      fiber_ind,
                      potentials_path,
                      waveform_path,
            )


def make_run_sub_list(run_number: int):
    """Create a list of all fiber simulations to be run. Skips fiber sims with existing output.
    :param run_number: the number of the run
    :return: a dict of all fiber simulations to be run
    """
    # build configuration filename
    filename: str = os.path.join('runs', f'{run_number}.json')
    # load in configuration data
    run = load(filename)

    submit_list = {}

    # assign appropriate configuration data
    samples = [run.get('sample', [])]
    models = run.get('models', [])
    sims = run.get('sims', [])

    for sample in samples:
        # loop models, sims
        for model in models:
            for sim in sims:
                sim_dir = os.path.join('n_sims')
                sim_name_base = '{}_{}_{}_'.format(sample, model, sim)
                nsim_list = [x for x in os.listdir(sim_dir) if x.startswith(sim_name_base)]
                for sim_name in nsim_list:
                    submit_list[sim_name] = []

                    sim_path = os.path.join(sim_dir, sim_name)
                    fibers_path = os.path.abspath(os.path.join(sim_path, 'data', 'inputs'))
                    output_path = os.path.abspath(os.path.join(sim_path, 'data', 'outputs'))

                    n_sim = sim_name.split('_')[-1]
                    sim_config = load(os.path.join(sim_path, '{}.json'.format(n_sim)))

                    fibers_files = [x for x in os.listdir(fibers_path) if re.match('inner[0-9]+_fiber[0-9]+\\.dat', x)]

                    for i, fiber_filename in enumerate(fibers_files):
                        master_fiber_name = str(fiber_filename.split('.')[0])
                        inner_name, fiber_name = tuple(master_fiber_name.split('_'))
                        inner_ind = int(inner_name.split('inner')[-1])
                        fiber_ind = int(fiber_name.split('fiber')[-1])

                        if sim_config['protocol']['mode'] == 'FINITE_AMPLITUDES':
                            n_amp = len(sim_config['protocol']['amplitudes'])
                            search_path = os.path.join(
                                output_path,
                                'activation_inner{}_fiber{}_amp{}.dat'.format(inner_ind, fiber_ind, n_amp - 1),
                            )
                        else:
                            search_path = os.path.join(
                                output_path,
                                f"thresh_inner{inner_ind}_fiber{fiber_ind}.dat",
                            )

                        if os.path.exists(search_path):
                            if args.verbose:
                                print(
                                    'Found {} -->\t\tskipping inner ({}) fiber ({})'.format(
                                        search_path, inner_ind, fiber_ind
                                    )
                                )
                                time.sleep(1)
                            continue

                        submit_list[sim_name].append({"job_number": i, "inner": inner_ind, "fiber": fiber_ind})
                    # save_submit list as csv
                    pd.DataFrame(submit_list[sim_name]).to_csv(os.path.join(sim_path, 'out_err_key.csv'), index=False)

    return submit_list

def confirm_submission(n_fibers, rundata, submission_context):
    """Confirm that the user wants to submit the simulations.
    :param n_fibers: the number of fibers to be run
    :param rundata: the run data (JSON config)
    :param submission_context: the submission context (e.g. cluster or local)
    """
    if n_fibers == 0:
        sys.exit('No fibers to run. Exiting...')
    if not args.skip_summary:
        # format run data
        df = pd.DataFrame(rundata)
        df.RUN = df.RUN.astype(int)
        df = df.sort_values('RUN')
        # print out and check that the user is happy
        print('Submitting the following runs (submission_context={}):'.format(submission_context))
        print(df.to_string(index=False))
        print('Will result in running {} fiber simulations'.format(n_fibers))
        proceed = input('\t Would you like to proceed?\n' '\t\t 0 = NO\n' '\t\t 1 = YES\n')
        if int(proceed) != 1:
            sys.exit()
        else:
            print('Proceeding...\n')
    else:
        print('Skipping summary, submitting {} fibers...'.format(n_fibers))


def get_submission_list(run_inds):
    """Get the list of simulations to be submitte for all runs.
    :param run_inds: the list of run indices
    :return: summary of runs, a list of all simulations to be submitted
    """
    rundata = []
    submission_list = {}
    for run_number in run_inds:

        # build configuration filename
        filename = os.path.join('runs', f'{run_number}.json')

        # configuration file exists
        assert os.path.exists(filename), 'Run configuration not found: {}'.format(run_number)

        # load in configuration data
        run = load(filename)

        # configuration is not empty
        assert len(run.items()) > 0, 'Encountered empty run configuration: {}'.format(filename)

        print('Generating run list for run {}'.format(run_number))
        # sleep to make it not too fast
        time.sleep(1)
        # get list of fibers to run
        submission_addition = make_run_sub_list(run_number)
        # check for duplicate nsims
        if any([x in submission_list for x in submission_addition.keys()]):
            warnings.warn('Duplicate nsims found in run {}. Continuing'.format(run_number))
        submission_list.update(submission_addition)
        rundata.append(
            {
                'RUN': run_number,
                'SAMPLE': run['sample'],
                'MODELS': run['models'],
                'SIMS': run['sims'],
            }
        )
    return rundata, submission_list


def pre_submit_setup():
    """Setup for submitting simulations.
    :return: the list of runs to be submitted, submission_context
    """
    # validate inputs
    global args
    args = parser.parse_args()
    if args.all_runs is True:
        if len(args.run_indices) > 0:
            sys.exit('Error: Cannot use -A/--run-all argument and pass run indices.')
        args.run_indices = [int(os.path.splitext(file)[0]) for file in os.listdir('runs') if file.endswith('.json')]
    if len(args.run_indices) == 0:
        sys.exit("Error: No run indices to use.")
    run_inds = args.run_indices
    # compile MOD files if they have not yet been compiled
    auto_compile(args.force_recompile)
    # check for submission context
    if args.cluster_submit:
        submission_context = 'cluster'
    elif args.local_submit:
        submission_context = 'local'
    else:
        submission_context = 'cluster' if shutil.which('sbatch') is not None else 'local'

    return run_inds, submission_context


# main
def main():
    """Main function."""
    # pre submit setup
    run_inds, submission_context = pre_submit_setup()
    # get list of simulations to be submitted
    rundata, submission_list = get_submission_list(run_inds)
    # confirm that the user wants to submit the simulations
    n_fibers = sum([len(x) for x in submission_list.values()])
    confirm_submission(n_fibers, rundata, submission_context)
    # make shell scripts for fiber submission
    make_fiber_tasks(submission_list, submission_context)
    # submit fibers
    submit_fibers(submission_context, submission_list)

if __name__ == "__main__":  # Allows for the safe importing of the main module
    main()
    print('done')
