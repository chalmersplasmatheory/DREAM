#!/usr/bin/python3 -i

import argparse
from DREAM import DREAMOutput, setup_interactive, who
import numpy as np
import os
import sys


def collect_single(name, output):
    """
    Collect data for the unknown quantity named 'name' from the list
    of DREAMOutput objects 'output'.
    """
    data = []
    for do in output:
        v = do.eqsys[name].data[0,:]

        data.append(v)

    return np.array(data)


def load_quantities(directory, timestep):
    """
    Load debug output from the named directory and assemble into a single
    DREAMOutput object.
    
    :param pattern: String with the format
    """
    pattern = f'{directory}/debugout_{{}}_{{}}.h5'

    do = DREAMOutput(pattern.format(timestep, 1))

    # Load all output objects
    dos = [do]
    iteration = 1
    while os.path.exists(pattern.format(timestep, iteration)):
        dos.append(DREAMOutput(pattern.format(timestep, iteration)))
        iteration += 1

    # Extend main output object
    do.grid.t = np.array(range(len(dos))) + 1
    for variable, uqn in do.eqsys.unknowns.items():
        data = collect_single(variable, dos)

        attr = {'description': uqn.description, 'equation': uqn.description_eqn}
        do.eqsys.setUnknown(name=variable, data=data, attr=attr)

    # Close inputs (except 'do', which is index 0)
    for i in range(1, len(dos)):
        dos[i].close()

    return do


def create_argparser():
    parser = argparse.ArgumentParser(description="DREAM Debug Output CLI")

    parser.add_argument('-t', '--timestep', help="Timestep index to load", type=int, default=0)
    parser.add_argument('directory', help="Name of directory containing debug output files to load", type=str, default=".")

    return parser.parse_args()


def determine_timestep(directory):
    """
    Determine which time step to load from the named directory.
    """
    pattern = f'{directory}/debugout_{{}}_{{}}.h5'
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

    for f in os.listdir(directory):
        if f.startswith('debugout_') and f.endswith('.h5'):
            return int(f.split('_')[1])

    raise Exception("Failed to automatically determine which time step to load.")
    

def main():
    args = create_argparser()

    if args.timestep <= 0:
        timestep = determine_timestep(args.directory)
        print(f'Loading time step {timestep}')
    else:
        timestep = args.timestep

    do = load_quantities(args.directory, timestep)
    setup_interactive(do, glob=globals())


if __name__ == '__main__':
    main()


