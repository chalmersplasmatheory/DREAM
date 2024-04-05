import glob
import matplotlib.pyplot as plt

import DREAM



def plot_charge_state_densities(ax, do):
    cmap = plt.cm.get_cmap("plasma_r", do.eqsys.n_i.getMultiples())
    nre = do.eqsys.n_re.data[-1,:]
    for i, ion in enumerate(do.eqsys.n_i.ions[0]):
        ax.semilogx(nre, ion.data[-1,:]/1e22, c=cmap(i), label=ion.name)
    # ax.legend()

def plot_electron_density(ax, do):
    nre = do.eqsys.n_re.data[-1,:]
    ncold = do.eqsys.n_cold.data[-1,:]
    ax.semilogx(nre, ncold)

def plot_effective_charge(ax, do):
    nre = do.eqsys.n_re.data[-1,:]
    Zeff = do.other.fluid.Zeff.data[-1,:]
    ax.semilogx(nre, Zeff)




if __name__ == '__main__':


    nre, ncold, Zeff = list(), list(), list()
    for file in glob.glob("outputs/kinetic*.h5"):

        do = DREAM.DREAMOutput(file)

        nre.append(do.eqsys.n_re.data[-1,0])
        ncold.append(do.eqsys.n_cold.data[-1,0])
        Zeff.append(do.other.fluid.Zeff.data[-1,0])

    fig, (ax1, ax2) = plt.subplots(ncols=2)

    ax1.scatter(nre, ncold)
    ax2.scatter(nre, Zeff)

    do = DREAM.DREAMOutput("outputs/fluid.h5")
    plot_electron_density(ax1, do)
    plot_effective_charge(ax2, do)

    plt.show()
