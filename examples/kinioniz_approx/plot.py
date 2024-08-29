#!/usr/bin/env python3
import glob
import matplotlib.pyplot as plt

import DREAM.DREAMOutput

if __name__ == '__main__':


    nre, ncold, Zeff = list(), list(), list()
    for file in glob.glob("outputs/kinetic*.h5"):

        do = DREAM.DREAMOutput(file)

        nre.append(do.eqsys.n_re.data[-1,0])
        ncold.append(do.eqsys.n_cold.data[-1,0])
        Zeff.append(do.other.fluid.Zeff.data[-1,0])

    fig, (ax1, ax2) = plt.subplots(ncols=2)

    ax1.scatter(nre, ncold, label="kinetic")
    ax2.scatter(nre, Zeff)

    do = DREAM.DREAMOutput("outputs/fluid.h5")
    ax1.semilogx(do.eqsys.n_re.data[-1,:], do.eqsys.n_cold.data[-1,:], label="fluid")
    ax2.semilogx(do.eqsys.n_re.data[-1,:], do.other.fluid.Zeff.data[-1,:])

    ax1.legend()
    ax1.set_ylabel("electron density [m$^{-3}$]")
    ax1.set_xlabel("runaway density [m$^{-3}$]")
    ax2.set_ylabel("effective charge [e]")
    ax2.set_xlabel("runaway density [m$^{-3}$]")
    

    plt.tight_layout()
    plt.show()
