import os
import numpy as np
import matplotlib.pyplot as plt

import DREAM.Settings.Equations.IonSpecies as Ions
import DREAM.Settings.Equations.ColdElectronTemperature as Temperature
from DREAM import DREAMSettings, DREAMOutput, runiface

N_POINTS = 1_000
ELECTRON_TEMPERATURE_LOW = 1e0
ELECTRON_TEMPERATURE_HIGH = 1e5
ELECTRON_DENSITY = 1e20



def setup_common():
    
    # DREAM simulation setup
    ds = DREAMSettings()
    ds.other.include(["fluid", "scalar"])
    
    # no kinetic grids needed
    ds.runawaygrid.setEnabled(False)
    ds.hottailgrid.setEnabled(False)

    # set geometry
    ds.radialgrid.setB0(5)
    ds.radialgrid.setMinorRadius(1)
    ds.radialgrid.setWallRadius(1.1)
    ds.radialgrid.setNr(N_POINTS)

    # set electric field
    ds.eqsys.E_field.setPrescribedData(1e-4)

    # set electron temperature
    r = np.linspace(0, 1, N_POINTS+1)[1:]
    Te = np.logspace(np.log10(ELECTRON_TEMPERATURE_LOW), np.log10(ELECTRON_TEMPERATURE_HIGH), N_POINTS)
    ds.eqsys.T_cold.setInitialProfile(Te, radius=r)
    ds.eqsys.T_cold.setType(Temperature.TYPE_SELFCONSISTENT) 
    
    return ds


def run_simulation(cW, include_recomb=True):
    """
    Sets up a simple DREAM simulation with or without recombination radiation.
    Returns the corresponding DREAMOutput object.
    
    float cW:               Tungsten concentration cW = nW/ne.
    bool include_recomb:    If 'False', neglects any recombination radiation.
    """
    filename = f"output_recomb={include_recomb}.h5"
    if os.path.exists(filename):
        #return DREAMOutput(filename)
        os.remove(filename)
        os.remove("test.h5")

    ## init: obtain coronal equilibrium distributions (fix ne and cW)

    ds = setup_common()

    # set deuterium and tungsten ions
    ds.eqsys.n_i.addIon(name="D", Z=1, isotope=2, init_equil=True, iontype=Ions.IONS_DYNAMIC, n=ELECTRON_DENSITY)
    if cW > 0:
        ds.eqsys.n_i.addIon(name="W", Z=74, init_equil=True, iontype=Ions.IONS_DYNAMIC, n=cW*ELECTRON_DENSITY)

    # run init simulation
    ds.timestep.setTmax(1e-10)
    ds.timestep.setNt(1)
    do_init = runiface(ds, outfile="test.h5")
    if cW == 0:
        return do_init

    # adjust charge states
    nD = do_init.eqsys.n_i.getIonByName("D")
    nW = do_init.eqsys.n_i.getIonByName("W")
    phiD1 = nD.ionstates[1].data[-1,:] / nD.data[-1,...].sum(axis=0)
    nD_new = ELECTRON_DENSITY / phiD1
    for nWj in nW.ionstates[1:]:
        phiWj = nWj.data[-1,:] / nW.data[-1,...].sum(axis=0)
        nD_new -= ELECTRON_DENSITY * cW * nWj.Z0 * phiWj / phiD1
   

    ## main: run with correct ion setup
    
    ds = setup_common()
    
    if include_recomb:
        ds.eqsys.T_cold.setRecombinationRadiation(Temperature.RECOMBINATION_RADIATION_INCLUDED)
    else:
        ds.eqsys.T_cold.setRecombinationRadiation(Temperature.RECOMBINATION_RADIATION_NEGLECTED)

    # set deuterium and tungsten ions
    ds.eqsys.n_i.addIon(name="D", Z=1, isotope=2, iontype=Ions.IONS_DYNAMIC, Z0=1, n=nD_new, r=do_init.grid.r)
    ds.eqsys.n_i.addIon(name="W", Z=74, init_equil=True, iontype=Ions.IONS_DYNAMIC, n=cW*ELECTRON_DENSITY)

    # run main simulation
    ds.timestep.setTmax(1e-10)
    ds.timestep.setNt(1)
    do = runiface(ds, outfile=filename)

    return do



def plot_Prad(do, ax):
    """
    Plots radiated power densities as function of temperature.

    DREAMOutput do: DREAM simulation output object.
    plt.Axes ax:    Pyplot Axes object used for plotting.
    """
    Te = do.eqsys.T_cold.data[-1,:]
    Prad_tot = do.other.fluid.Tcold_radiation.data[-1,:]
    Prad_D = do.other.fluid.radiated_power_ion.data[-1,0,:]
    Prad_W = do.other.fluid.radiated_power_ion.data[-1,1,:]

    ax.plot(Te, Prad_tot, 'k--', label="total")
    ax.plot(Te, Prad_D, 'r', label="deuterium")
    ax.plot(Te, Prad_W, 'b', label="tungsten")


def plot_Pion(do, ax):
    """
    Plots change in potential energy due to ionisation and recombination as function of temperature.

    DREAMOutput do: DREAM simulation output object.
    plt.Axes ax:    Pyplot Axes object used for plotting.
    """
    Te = do.eqsys.T_cold.data[-1,:]
    Pion_tot = do.other.fluid.Tcold_binding_energy.data[-1,:]
    Pion_D = do.other.fluid.binding_energy_ion.data[-1,0,:]
    Pion_W = do.other.fluid.binding_energy_ion.data[-1,1,:]

    ax.plot(Te, Pion_tot, 'k--', label="total")
    ax.plot(Te, Pion_D, 'r', label="deuterium")
    ax.plot(Te, Pion_W, 'b', label="tungsten")





def main():


    import argparse

    parser = argparse.ArgumentParser(description="Compare Prad and Pion with the contributions from individual ions (D and W).")
    parser.add_argument("--include_recombination", action="store_true", help="Include recombination radiation")
    parser.add_argument("--tungsten_concentration", type=float, default=1e-6, help="Tungsten concentration nW/ne")
    args = parser.parse_args()

    assert 0 <= args.tungsten_concentration < 1


    do = run_simulation(args.tungsten_concentration, include_recomb=args.include_recombination)


    


    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

    plot_Prad(do, ax1)
    plot_Pion(do, ax2)

    
    ax1.set_ylabel(r"$P_{\rm rad}$ [W/m$^3$]")
    ax2.set_ylabel(r"$P_{\rm ion}$ [W/m$^3$]")
    ax2.set_xlabel(r"$T_e$ [eV]")
    ax1.set_yscale("log")
    ax2.set_yscale("symlog")
    ax2.set_xscale("log")
    ax1.legend()
    
    plt.subplots_adjust(hspace=0)
    plt.show()

if __name__ == "__main__":
    main()


