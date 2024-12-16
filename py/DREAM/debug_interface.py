#!/usr/bin/python3 -i

import argparse
from DREAM import DREAMOutput, setup_interactive, who
from DREAM.Output import UnknownQuantity
import matplotlib.pyplot as plt
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


def load_quantities(directory, timestep, outputobject=DREAMOutput):
    """
    Load debug output from the named directory and assemble into a single
    DREAMOutput object.
    
    :param pattern: String with the format
    """
    pattern = f'{directory}/debugout_{{}}_{{}}.h5'

    do = outputobject(pattern.format(timestep, 1), loadsettings=True)

    # Load all output objects
    dos = [do]
    iteration = 1
    while os.path.exists(pattern.format(timestep, iteration)):
        dos.append(outputobject(pattern.format(timestep, iteration), loadsettings=False))
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


def plotConvergenceOutput(do):
    """
    Plot the convergence of the unknowns in the given DREAMOutput object.
    """
    fig, ax = plt.subplots(1,1)

    lines = []
    for uname, u in do.eqsys.unknowns.items():
        reltol = do.settings.solver.tolerance.getRelTol(uname)

        if np.all(u.data==0):
            continue

        # Delta = np.abs((u.data[:]/ u.data[-1,:])-1)
        temp = np.sqrt(np.sum(u.data[:]**2,axis=tuple(range(1,u.data.ndim))))
        Delta = temp/(temp[-1])-1

        v = np.sum(Delta) / retol
        l, = ax.semilogy(do.grid.t, v, label=uname)
        lines.append(l)

    ax.plot([1, do.grid.t.size], [1, 1], 'k--')
    ax.set_xlim([1, do.grid.t.size])

    annot = ax.annotate("", xy=(0,0), xytext=(-20,20), textcoords="offset points", bbox={'boxstyle': 'round', 'fc':'w'}, arrowprops={'arrowstyle':'->'})
    annot.set_visible(False)

    def update_annot(line, ind):
        x, y = line.get_data()
        annot.xy = (x[ind['ind'][0]], y[ind['ind'][0]])
        #annot.set(backgroundcolor=line.get_color())
        annot.set(bbox={'boxstyle':'round', 'fc': 'w', 'ec': line.get_color()})
        annot.set_text(line.get_label())
    
    def hover(event):
        if event.inaxes == ax:
            found = False
            for line in lines:
                cont, ind = line.contains(event)

                if cont:
                    update_annot(line, ind)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                    found = True
                    break

            if not found and annot.get_visible():
                annot.set_visible(False)
                fig.canvas.draw_idle()

    fig.canvas.mpl_connect('motion_notify_event', hover)

    plt.show(block=False)


def plotConvergenceUnknown(u):
    """
    Plot the convergence of the given unknown quantity.
    """
    pass


def plotConvergence(u=None):
    """
    Plot the convergence of a given unknown quantity of DREAMOutput object.
    """
    if u == None:
        u = glob['do']
    
    if isinstance(u, type(UnknownQuantity)):
        plotConvergenceUnknown(u)
    else:
        plotConvergenceOutput(u)


def create_argparser():
    parser = argparse.ArgumentParser(description="DREAM Debug Output CLI")

    parser.add_argument('-t', '--timestep', help="Timestep index to load", type=int, default=0)
    parser.add_argument('directory', help="Name of directory containing debug output files to load", nargs='?', type=str, default=".")

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
    

def main(glob=None, outputobject=DREAMOutput):
    if glob is None:
        glob = globals()

    args = create_argparser()

    if args.timestep <= 0:
        timestep = determine_timestep(args.directory)
        print(f'Loading time step {timestep}')
    else:
        timestep = args.timestep

    do = load_quantities(args.directory, timestep, outputobject=outputobject)
    globals()['glob'] = glob
    setup_interactive(do, glob=glob)


if __name__ == '__main__':
    main(glob=globals())


