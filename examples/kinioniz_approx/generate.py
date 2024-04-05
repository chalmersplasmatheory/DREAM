import numpy as np
import matplotlib.pyplot as plt
import scipy.constants

import DREAM
import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.DistributionFunction as DistFunc
import DREAM.Settings.Equations.RunawayElectronDistribution as REDist

def run(ds, tmax, file="output.h5"):
    ds.timestep.setIonization(dt0=1e-7, dtmax=1e-6, tmax=tmax)
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

    ds.eqsys.n_i.addIon(name='Ar', Z=18, n=1e22, iontype=Ions.IONS_DYNAMIC_FULLY_IONIZED)


    ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_FLUID_APPROX_RE)


    ds.other.include(["fluid/Zeff"])

    return ds


def getAvalancheDistribution(p, xi, E, Z, nre=1, logLambda=15):
    """
    Evaluates the analytical avalanche distribution function according to
    equation (2.17) of [Embreus et al, J. Plasma Phys. 84 (2018)].

    :param p:         Momentum grid on which to evaluate the distribution.
    :param xi:        Pitch grid on which to evaluate the distribution.
    :param E:         Electric field strength (normalized to the Connor-Hastie field, Ec).
    :param Z:         Plasma total charge (= 1/ne_tot * sum_i ni * Zi^2)
    :param nre:       Runaway electron density.
    :param logLambda: Coulomb logarithm.
    """
    if p.ndim == 1:
        P, XI = np.meshgrid(p, xi)
    else:
        P, XI = p, xi

    c = scipy.constants.c
    m_e = scipy.constants.m_e

    g = np.sqrt(1+P**2)
    A = (E+1) / (Z+1) * g
    cZ = np.sqrt(5+Z)
    g0 = cZ*logLambda

    pf = m_e*c * nre * A / (2*np.pi*m_e*c*g0*P**2) / (1-np.exp(-2*A))
    f = pf * np.exp(-g/g0 - A*(1-XI))

    return f



def generate_kinetic(nre, Tcold=1e4):
    ds = generate_fluid(nre, Tcold)
    ds.eqsys.n_i.setIonization(Ions.IONIZATION_MODE_KINETIC)

    Np = 100
    Nxi = 40
    pMin, pMax = 1, 100

    ds.runawaygrid.setEnabled(True)
    ds.runawaygrid.setNp(Np)
    ds.runawaygrid.setNxi(Nxi)
    ds.runawaygrid.setPmin(pMin)
    ds.runawaygrid.setPmax(pMax)

    # Distribution function
    f   = np.zeros((1, ds.radialgrid.nr, Nxi, Np))
    pp  = np.linspace(pMin, pMax, Np+1)
    xip = np.linspace(-1, 1, Nxi+1)
    p   = 0.5 * (pp[1:] + pp[:-1])
    xi  = 0.5 * (xip[1:] + xip[:-1])

    # Avalanche distribution parameters
    Ztot = 8 # Total plasma charge (c.f. [Embreus JPP 84 (2018)])


    for ir in range(len(nre)):
        f[0,ir,:] = getAvalancheDistribution(p=p, xi=xi, E=1, Z=Ztot, nre=nre[ir])

    r = np.linspace(0, 1, len(nre))
    ds.eqsys.f_re.prescribe(f=f, t=[0], r=r, xi=xi, p=p)

    ds.eqsys.f_re.setAdvectionInterpolationMethod(ad_int=DistFunc.AD_INTERP_TCDF)



    return ds


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
    # ds = generate_fluid(np.logspace(14, 21.27, 400), Tcold=1)
    # do = run(ds, 1e-2, file="fluid.h5")

    ds = generate_kinetic(np.logspace(14, 19, 20), Tcold=1)
    do = run(ds, 1e-2, file="kinetic.h5")

    # do = DREAM.DREAMOutput("kinetic.h5")

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)
    plot_charge_state_densities(ax1, do)
    plot_electron_density(ax2, do)
    plot_effective_charge(ax3, do)
    plt.show()
