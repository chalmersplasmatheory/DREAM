/**
 * A test for the BootstrapCurrent term in DREAM.
 */

#include <vector>
#include <string>
#include <gsl/gsl_interp.h>
#include "BootstrapCurrent.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/Equations/BootstrapCurrent.hpp"
#include "DREAM/Equations/Fluid/BootstrapElectronDensityTerm.hpp"
#include "DREAM/Equations/Fluid/BootstrapElectronTemperatureTerm.hpp"
#include "DREAM/Equations/Fluid/BootstrapIonDensityTerm.hpp"
#include "DREAM/Equations/Fluid/BootstrapIonThermalEnergyTerm.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

#include <iostream>

using namespace DREAMTESTS::_DREAM;
using namespace std;

const string INPUT_FILENAME = "/home/ida/Desktop/DREAM/tests/cxx/tests/DREAM/profiles_IDE_40655_t=2.3s.h5";
const string EQUIL_FILENAME = "/home/ida/Desktop/DREAM/tests/cxx/tests/DREAM/equilibrium_IDE_40655_t=2.3s.h5";

const len_t N_IONS = 2; // Number of ion species
const len_t Z_IONS[N_IONS] = {1, 10}; // Vector with atomic number for each ion species, # elements == N_IONS
const char ION_NAMES[N_IONS][3] = {"D", "Ne"}; // Vector with name of each ion species, # elements == N_IONS


/**
 * Generate a default ion handler.
 */
DREAM::IonHandler *BootstrapCurrent::GetIonHandler(
    DREAM::FVM::Grid *g, DREAM::FVM::UnknownQuantityHandler *uqh
) {
    vector<string> tritiumNames(0), hydrogenNames(0); 
    vector<string> names(N_IONS);
    len_t *Z = new len_t[N_IONS];// The ion charge numbers must be provided to the IONHandler as a dynamically allocated array to avoid memory issues
    for (len_t i = 0; i < N_IONS; i++){
        names[i] = ION_NAMES[i];
        Z[i]=Z_IONS[i];
    }

    return new DREAM::IonHandler(
        g->GetRadialGrid(), uqh, Z, N_IONS, names, tritiumNames, hydrogenNames
    );
}

/**
 * Generate an unknown quantity handler.
 */
DREAM::FVM::UnknownQuantityHandler *BootstrapCurrent::GetUnknownHandler(
	DREAM::FVM::Grid *g, bool withIonEnergy, SFile_HDF5 *input_file
) {
    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();

	len_t nr = g->GetNr();

    len_t nZ0 = 0;
    for (len_t i = 0; i < N_IONS; i++)
        nZ0 += Z_IONS[i] + 1;

    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_ION_SPECIES, "0", g, nZ0);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_COLD, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_T_COLD, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_NI_DENS, "0", g, N_IONS);
    if (withIonEnergy)
        uqh->InsertUnknown(DREAM::OptionConstants::UQTY_WI_ENER, "0", g, N_IONS);

    
    const real_t *r = g->GetRadialGrid()->GetR();
    gsl_interp_accel *acc_r = gsl_interp_accel_alloc();
	sfilesize_t dims[2];
    const real_t *input_r = input_file->GetDoubles1D("radius", dims);
    // Construct splines for input data
    // TODO Peter: How does the data look? Is it 1D, only radius? Or is it nt*nr? Change dims[0] accordingly
    gsl_spline *spline_ne  = gsl_spline_alloc(gsl_interp_steffen, dims[0]);
    gsl_spline_init(spline_ne, input_r, input_file->GetDoubles1D("electron_density", dims), dims[0]);

    gsl_spline *spline_Te  = gsl_spline_alloc(gsl_interp_steffen, dims[0]);
    gsl_spline_init(spline_Te, input_r, input_file->GetDoubles1D("electron_temperature", dims), dims[0]);

    gsl_spline *spline_ni  = gsl_spline_alloc(gsl_interp_steffen, dims[0]);
    gsl_spline_init(spline_ni, input_r, input_file->GetDoubles1D("ion_density", dims), dims[0]);

    gsl_spline *spline_Ti  = gsl_spline_alloc(gsl_interp_steffen, dims[0]);
    gsl_spline_init(spline_Ti, input_r, input_file->GetDoubles1D("ion_temperature", dims), dims[0]);

    gsl_spline *spline_Zeff  = gsl_spline_alloc(gsl_interp_steffen, dims[0]);
    gsl_spline_init(spline_Zeff, input_r, input_file->GetDoubles1D("effective_charge", dims), dims[0]);
	

    // Set initial values
    const len_t Nsum = N_IONS*g->GetNr();
    const len_t N = nZ0*g->GetNr();;
    real_t *Nions = new real_t[Nsum];
    real_t *nions = new real_t[N];
    len_t  rOffset = 0;
	len_t Z = Z_IONS[1];
    for (len_t iIon = 0; iIon < N_IONS; iIon++) {
        for (len_t Z0 = 0; Z0 <= Z_IONS[iIon]; Z0++) {
            for (len_t ir = 0; ir < nr; ir++, rOffset++){
                real_t ni = gsl_spline_eval(spline_ni, r[ir], acc_r);
                real_t Zeff = gsl_spline_eval(spline_Zeff, r[ir], acc_r);
				if (Z0 == 1 && iIon == 0) // deuterium
					nions[rOffset] = ni * Z*(Z - Zeff) / ((Z - 1)*(Z - Zeff + 1));
				else if (Z0 == Z_IONS[1] && iIon == 1) // impurity
					nions[rOffset] = ni * (Zeff - 1) / (Z*Z - Zeff*(Z - 1) - 1);
				else
					nions[rOffset] = 0.0;
				Nions[iIon*nr+ir] += nions[rOffset]; // only one charge state per ion
			}
		}
	}
            
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_ION_SPECIES, nions);
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_NI_DENS, Nions);

    if (withIonEnergy){
        real_t *wions = new real_t[Nsum];
        len_t rOffset = 0;
        for (len_t iIon = 0; iIon < N_IONS; iIon++) {
            for (len_t ir = 0; ir < nr; ir++, rOffset++) {
                real_t Ti = gsl_spline_eval(spline_Ti, r[ir], acc_r);
                wions[rOffset] = Ti * 1.5 * DREAM::Constants::ec * Nions[iIon*nr+ir];
            }
        }
        uqh->SetInitialValue(DREAM::OptionConstants::UQTY_WI_ENER, wions);

        delete [] wions;
    }

    // Set electron quantities
    real_t *ne = new real_t[g->GetNr()];
    real_t *Te = new real_t[g->GetNr()];
    for (len_t ir = 0; ir < nr; ir++) {
        ne[ir] = gsl_spline_eval(spline_ne, r[ir], acc_r);
        Te[ir] = gsl_spline_eval(spline_Te, r[ir], acc_r);
    }
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_N_COLD, ne);
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_T_COLD, Te);
    
    gsl_spline_free(spline_ne);
    gsl_spline_free(spline_Te);
    gsl_spline_free(spline_ni);
    gsl_spline_free(spline_Ti);
    gsl_spline_free(spline_Zeff);
    gsl_interp_accel_free(acc_r);

    delete [] nions;
    delete [] Nions;
    delete [] ne;
    delete [] Te;
    
    return uqh;
}

/**
 * Check if the bootstrap current is consistent with IDA.
 */
bool BootstrapCurrent::CheckBootstrap(bool withIonEnergy) {

	// load input data
	SFile_HDF5 *sf = new SFile_HDF5();
    
    sf->Open(INPUT_FILENAME, SFILE_MODE_READ);

	sfilesize_t dims[2];
	real_t *j_bs_IDA_f = sf->GetDoubles1D("bootstrap_current", dims);

	
    real_t successRelErrorThreshold = 1e-5; // TODO Peter: Probably want to change this
    bool success = true;
    const len_t Nr = (len_t)dims[0]-1;	// match radial grid with that of the input data 
    DREAM::FVM::Grid *grid = this->InitializeNumericFluidGrid(sf->GetDoubles1D("radius", dims), Nr, EQUIL_FILENAME);
    DREAM::FVM::UnknownQuantityHandler *uqh = GetUnknownHandler(grid, withIonEnergy, sf);
    DREAM::IonHandler *ih = GetIonHandler(grid, uqh);
    ih->Rebuild();

    DREAM::OptionConstants::momentumgrid_type gridtype = DREAM::OptionConstants::MOMENTUMGRID_TYPE_PXI;
    DREAM::CollisionQuantity::collqty_settings *cq = new DREAM::CollisionQuantity::collqty_settings;

    cq->collfreq_type = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_PARTIALLY_SCREENED;
    cq->collfreq_mode = DREAM::OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
    cq->lnL_type      = DREAM::OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;
    cq->bremsstrahlung_mode = DREAM::OptionConstants::EQTERM_BREMSSTRAHLUNG_MODE_STOPPING_POWER;

    DREAM::CoulombLogarithm *lnLEE = new DREAM::CoulombLogarithm(grid,uqh,ih,gridtype,cq,DREAM::CollisionQuantity::LNLAMBDATYPE_EE);
    DREAM::BootstrapCurrent *bootstrap = new DREAM::BootstrapCurrent(grid, uqh, ih, lnLEE);
    bootstrap->Rebuild();

    // Construct solution vector
    real_t *j_bs = new real_t[Nr];
    real_t *j_bs_IDA = new real_t[Nr];
    
    // Clear 'j_bs'
    for (len_t ir = 0; ir < Nr; ir++){
        j_bs[ir] = 0;
        j_bs_IDA[ir] = (j_bs_IDA_f[ir+1]+j_bs_IDA_f[ir])/2;
    }
    
    // Construct equation for each ion species
    DREAM::BootstrapElectronDensityTerm *term_ncold = new DREAM::BootstrapElectronDensityTerm(grid, uqh, bootstrap, ih);
    DREAM::BootstrapElectronTemperatureTerm *term_Tcold = new DREAM::BootstrapElectronTemperatureTerm(grid, uqh, bootstrap, ih);

    DREAM::BootstrapIonDensityTerm *term_Ni[N_IONS];
    DREAM::BootstrapIonThermalEnergyTerm *term_wi[N_IONS];
    for (len_t iIon = 0; iIon < N_IONS; iIon++){
        term_Ni[iIon] = new DREAM::BootstrapIonDensityTerm(grid, uqh, bootstrap, ih, iIon);
        if (withIonEnergy)
            term_wi[iIon] = new DREAM::BootstrapIonThermalEnergyTerm(grid, uqh, bootstrap, ih, iIon);
    }

    
    term_ncold->Rebuild(0, 0, uqh);
    term_ncold->SetVectorElements(j_bs, nullptr);
    term_Tcold->Rebuild(0, 0, uqh);
    term_Tcold->SetVectorElements(j_bs, nullptr);
    for (len_t iIon = 0; iIon < N_IONS; iIon++) {
        term_Ni[iIon]->Rebuild(0, 0, uqh);
        term_Ni[iIon]->SetVectorElements(j_bs, nullptr);
        if (withIonEnergy){
            term_wi[iIon]->Rebuild(0, 0, uqh);
            term_wi[iIon]->SetVectorElements(j_bs, nullptr);
        }
    }

    // Sum of all elements should vanish
    real_t deltas;
    for (len_t ir = 0; ir < Nr; ir++) {
        deltas = abs(j_bs[ir] - j_bs_IDA[ir]) / j_bs_IDA[ir];
        printf("ir=%ld  --  DREAM: %.1f,       IDA: %.1f\n", ir, j_bs[ir], j_bs_IDA[ir]);
        if(deltas>successRelErrorThreshold){
            success = false;
            std::cout << "delta [rel error]: " << deltas << endl;
        }
    }


    // Deallocate equation terms
    delete term_ncold;
    delete term_Tcold;
    for (len_t iIon = 0; iIon < N_IONS; iIon++){
        delete term_Ni[iIon];
        delete term_wi[iIon];
    }
    delete [] j_bs;
    delete [] j_bs_IDA;
    delete bootstrap;
    delete lnLEE;
    delete cq;
    delete ih;
    delete uqh;
    delete grid;
    delete sf;

    return success;
}

/**
 * Run this test.
 */
bool BootstrapCurrent::Run(bool) {
    bool success = true;

    if ((success &= CheckBootstrap(true)) && (success &= CheckBootstrap(false)))
        this->PrintOK("The bootstrap current is consistent with IDA.");
    else
        this->PrintError("The bootstrap current is not consistent with IDA.");

    return success;
}

