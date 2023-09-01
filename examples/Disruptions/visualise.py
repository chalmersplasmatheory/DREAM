#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import DREAM


def load_outputs(filename_prefix):
    """
    Loads DREAMOutput objects from given load file pattern. If any is not found,
    it returns False, otherwise it returns a tuple with three DREAMOutput objects
    with filename suffixes '_init.h5', '_1.h5' and '_2.h5'.

    str filename_prefix:    Pattern to look for when loading the output files, eg.
                            'output/out_fluid_scen1_output_1'.
    """
    assert isinstance(filename_prefix, str), "expected a string as input"

    try:
        # do_init = DREAM.DREAMOutput(f"{filename_prefix}_init.h5")
        do_tq   = DREAM.DREAMOutput(f"{filename_prefix}_1.h5")
        do_cq   = DREAM.DREAMOutput(f"{filename_prefix}_2.h5")
    except FileNotFoundError:
        print(f"could not load from {filename_prefix}")
        return False

    return do_tq, do_cq


def concatenate_data(*dataobjects, current=False):
    """
    Concatenates data from given Dataobject objects into one single list without
    overlap in time, by removing the first element in time for each time point.
    Returns a tuple of two lists, one being time and the other the corresponding data.

    DREAM.Dataobject.Dataobject dataobject:    Data to concatenate.
    """
    assert len(dataobjects) > 1, "expected at least two inputs"
    assert all(obj.data.shape[1] == dataobjects[0].data.shape[1] for obj in dataobjects), "same no. radial points is required"

    time = dataobjects[0].time
    data = dataobjects[0].current() if current else dataobjects[0].data
    for i, obj in enumerate(dataobjects):
        if i == 0:
            continue
        time = np.concatenate((time, obj.time[1:] + time[-1]))
        data = np.concatenate((data, obj.current()[1:] if current else obj.data[1:,:]))

    # time = np.concatenate([obj.time[1:] for obj in dataobjects])
    # data = np.concatenate([obj.data[1:,:] for obj in dataobjects])

    return time * 1e6, data


def plot_temperature(filename_prefix, ax, c='k', lw=1):
    tT, T = concatenate_data(*[do.eqsys.T_cold for do in load_outputs(filename_prefix)])
    ax.semilogy(tT, T[:,4], c, lw=lw) # ir = 4 => r/a = .3

def plot_currents(filename_prefix, ax, c='k', lw=1):
    do_tq, do_cq = load_outputs(filename_prefix)
    tjtot, jtot = concatenate_data(do_tq.eqsys.j_tot, do_cq.eqsys.j_tot, current=True)
    _, jre = concatenate_data(do_tq.eqsys.j_re, do_cq.eqsys.j_re, current=True)
    ax.plot(tjtot, jtot * 1e-3, c, lw=lw)
    ax.plot(tjtot, jre * 1e-3, f"{c}--", lw=lw)


if __name__ == '__main__':

    fluid_sol_prefix = "output/out_fluid_sol_scen1_output"
    fluid_ure_prefix = "output/out_fluid_scen1_output"

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 6), sharex=True)

    plot_temperature(fluid_sol_prefix, ax1, c='b', lw=5)
    plot_currents(fluid_sol_prefix, ax2, c='b', lw=5)

    plot_temperature(fluid_ure_prefix, ax1, c='r', lw=1)
    plot_currents(fluid_ure_prefix, ax2, c='r', lw=1)

    plt.xlim(0, 300)
    ax1.set_ylabel("temperature [eV]")
    ax2.set_ylabel("current [kA]")
    ax2.set_xlabel(r"time [${\rm \mu s}$]")

    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    plt.show()
