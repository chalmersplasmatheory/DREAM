import numpy as np
import matplotlib.pyplot as plt

import DREAM
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.RunawayElectronDistribution as REDist

def run(ds, tmax, file="output.h5"):
    ds.timestep.setIonization(dt0=1e-7, dtmax=4e-6, tmax=tmax)
    ds.timestep.setMinSaveTimestep(3e-7)
    return DREAM.runiface(ds, outfile=file)


def generate_fluid(nre, Tcold=1e4):

    ds = DREAM.DREAMSettings()


    ds.hottailgrid.setEnabled(False)
    ds.runawaygrid.setEnabled(False)

    ds.radialgrid.setB0(1)
    ds.radialgrid.setMinorRadius(1)
    ds.radialgrid.setWallRadius(1.1)
    ds.radialgrid.setNr(10)

    if isinstance(nre, (int, float)):
        ds.eqsys.n_re.setInitialProfile(nre)
        ds.radialgrid.setNr(1)
    elif isinstance(nre, (list, np.ndarray)):
        ds.eqsys.n_re.setInitialProfile(nre, np.linspace(0, 1, len(nre)))
        ds.radialgrid.setNr(len(nre))


    ds.eqsys.E_field.setPrescribedData(1)
    ds.eqsys.T_cold.setPrescribedData(Tcold)

    ds.eqsys.n_i.addIon(name='Ar', Z=18, n=1e20, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED)


    ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_FLUID_APPROX_RE)


    ds.other.include(["fluid/Zeff"])

    return ds

def generate_kinetic(nre, Tcold=1e4):
    ds = generate_fluid(nre, Tcold)
    ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_KINETIC)

    ds.runawaygrid.setEnabled(True)
    ds.runawaygrid.setNxi(50)
    ds.runawaygrid.setNp(100)
    ds.runawaygrid.setPmax(30)
    # Use flux limiters
    ds.eqsys.f_re.setAdvectionInterpolationMethod(ad_int=DistFunc.AD_INTERP_TCDF)

    # Set initialization method
    ds.eqsys.f_re.setInitType(REDist.INIT_AVALANCHE)
    # ds.eqsys.f_re.enableAnalyticalDistribution(mode=DistFunc.DISTRIBUTION_MODE_PRESCRIBED)
    return ds


def plot_charge_state_densities(ax, do):
    cmap = plt.cm.get_cmap("plasma_r", do.eqsys.n_i.getMultiples())
    nre = do.eqsys.n_re.data[-1,:]
    for i, ion in enumerate(do.eqsys.n_i.ions[0]):
        ax.semilogx(nre, ion.data[-1,:]/1e20, c=cmap(i), label=ion.name)
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
    ds = generate_fluid(np.logspace(14, 21.27, 400), Tcold=1)
    do = run(ds, 1e-2, file="fluid.h5")

    ds = generate_kinetic(np.logspace(14, 21.27, 10), Tcold=1)
    do = run(ds, 1e-2, file="kinetic.h5")

    # do = DREAM.DREAMOutput("output.h5")

    # fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)
    # plot_charge_state_densities(ax1, do)
    # plot_electron_density(ax2, do)
    # plot_effective_charge(ax3, do)
    # plt.show()
