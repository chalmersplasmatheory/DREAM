/**
 * Implementation of ion equations.
 */

#include "DREAM/Equations/Fluid/IonPrescribedParameter.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/Equation.hpp"


using namespace DREAM;
using namespace std;


#define MODULENAME "equationsystem/n_i"

/**
 * Define options for the ions.
 *
 * set: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_Ions(Settings *s) {
    const len_t dims[1] = {0};

    s->DefineSetting(MODULENAME "/names", "Names of each ion species", (const string)"");
    s->DefineSetting(MODULENAME "/Z", "List of atomic charge numbers", 1, dims, (int_t*)nullptr);
    s->DefineSetting(MODULENAME "/types", "Method to use for determining ion charge distributions", 1, dims, (int_t*)nullptr);

    DefineDataIonR(MODULENAME, s, "initial");
    DefineDataIonRT(MODULENAME, s, "prescribed");
}

/**
 * Returns the number of ion charge states set by the configuration
 * (i.e. the number of elements "divided by nr (number of radial points)"
 * in the "ION_SPECIES" unknown quantity)
 *
 * set: Settings object to load charge state number from.
 */
len_t SimulationGenerator::GetNumberOfIonChargeStates(Settings *s) {
    len_t nZ;
    const int_t *Z = s->GetIntegerArray(MODULENAME "/Z", 1, &nZ, false);

    len_t nChargeStates = 0;
    for (len_t i = 0; i < nZ; i++)
        nChargeStates += Z[i]+1;
    
    return nChargeStates;
}

/**
 * Construct the equation governing the evolution of the
 * ion densities.
 */
void SimulationGenerator::ConstructEquation_Ions(EquationSystem *eqsys, Settings *s, ADAS *adas) {
    const real_t t0 = 0;
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    len_t nZ, ntypes;
    const string names = s->GetString(MODULENAME "/names");
    const int_t *_Z  = s->GetIntegerArray(MODULENAME "/Z", 1, &nZ);
    const int_t *itypes = s->GetIntegerArray(MODULENAME "/types", 1, &ntypes);

    // Parse list of ion names (stored as one contiguous string,
    // each substring separated by ';')
    vector<string> ionNames(std::count(names.begin(), names.end(), ';'));
    len_t si = 0;
    const len_t sl = names.size();
    for (len_t i = 0; i < sl; i++) {
        if (names[i] == ';') {
            si++;
            if (si == nZ && i+1 < sl)
                throw SettingsException(
                    "ions: Too many ion names given. Expected " LEN_T_PRINTF_FMT ".", nZ
                );
        } else
            ionNames[si] += names[i];
    }

    // Automatically name any unnamed ions
    if (ionNames.size() < nZ) {
        for (len_t i = ionNames.size(); i < nZ; i++)
            ionNames.push_back("Ion " + to_string(i));
    } else if (ionNames.size() > nZ) {
        throw SettingsException(
            "ions: Too many ion names given: %zu. Expected " LEN_T_PRINTF_FMT ".",
            ionNames.size(), nZ
        );
    }

    // Verify that exactly one type per ion species is given
    if (nZ != ntypes)
        throw SettingsException(
            "ions: Expected the lengths of 'Z' and 'types' to match."
        );

    // Data type conversion
    len_t *Z = new len_t[nZ];
    for (len_t i = 0; i < nZ; i++)
        Z[i] = (len_t)_Z[i];

    enum OptionConstants::ion_data_type *types = new enum OptionConstants::ion_data_type[ntypes];
    for (len_t i = 0; i < ntypes; i++)
        types[i] = (enum OptionConstants::ion_data_type)itypes[i];

    /////////////////////////
    /// LOAD ION DATA
    /////////////////////////
    // Count number of prescribed/dynamic charge states
    len_t nZ0_prescribed=0, nZ_prescribed=0, nZ_dynamic=0;
    len_t *prescribed_indices = new len_t[nZ];
    len_t *dynamic_indices = new len_t[nZ];
    for (len_t i = 0; i < nZ; i++) {
        switch (types[i]) {
            case OptionConstants::ION_DATA_PRESCRIBED:
                nZ0_prescribed += Z[i] + 1;
                prescribed_indices[nZ_prescribed++] = i;
                break;

            case OptionConstants::ION_DATA_TYPE_DYNAMIC:
            case OptionConstants::ION_DATA_EQUILIBRIUM:
                dynamic_indices[nZ_dynamic++] = i;
                break;

            default:
                throw SettingsException(
                    "ions: Unrecognized ion model type specified: %d.",
                    types[i]
                );
        }
    }

    // Load ion data
    real_t *dynamic_densities = LoadDataIonR(
        MODULENAME, fluidGrid->GetRadialGrid(), s, nZ_dynamic, "initial"
    );
    IonInterpolator1D *prescribed_densities = LoadDataIonRT(
        MODULENAME, fluidGrid->GetRadialGrid(), s, nZ0_prescribed, "prescribed"
    );

    IonHandler *ih = new IonHandler(fluidGrid->GetRadialGrid(), eqsys->GetUnknownHandler(), Z, nZ, ionNames);
    eqsys->SetIonHandler(ih);

    // Initialize ion equations
    FVM::Equation *eqn = new FVM::Equation(fluidGrid);

    IonPrescribedParameter *ipp = nullptr;
    if (nZ0_prescribed > 0)
        ipp = new IonPrescribedParameter(fluidGrid, ih, nZ_prescribed, prescribed_indices, prescribed_densities);

    // Construct dynamic equations
    for (len_t iZ = 0; iZ < nZ; iZ++) {
        switch (types[iZ]) {
            case OptionConstants::ION_DATA_PRESCRIBED: break;

            // 'Dynamic' and 'Equilibrium' differ by a transient term
            case OptionConstants::ION_DATA_TYPE_DYNAMIC:
                /*eqn->AddTerm(new IonTransientTerm(
                    fluidGrid, ih, iZ
                ));*/
                [[fallthrough]];

            case OptionConstants::ION_DATA_EQUILIBRIUM:
                eqn->AddTerm(new IonRateEquation(
                    fluidGrid, ih, iZ, adas, eqsys->GetUnknownHandler()
                ));
                break;

            default:
                throw SettingsException(
                    "ions: Unrecognized ion model type specified: %d.",
                    types[iZ]
                );
        }
    }

    if (ipp != nullptr)
        eqn->AddTerm(ipp);
    
    eqsys->SetEquation(OptionConstants::UQTY_ION_SPECIES, OptionConstants::UQTY_ION_SPECIES, eqn);

    // Initialize dynamic ions
    const len_t Nr = fluidGrid->GetNr();
    real_t *ni = new real_t[ih->GetNzs() * Nr];

    // Begin by evaluating prescribed densities
    if (ipp != nullptr)
        ipp->Evaluate(ni, nullptr);

    // ...and then fill in with the initial dynamic ion values
    for (len_t i = 0, ionOffset = 0; i < nZ_dynamic; i++) {
        len_t Z   = ih->GetZ(dynamic_indices[i]);
        len_t idx = ih->GetIndex(dynamic_indices[i], 0);

        for (len_t Z0 = 0; Z0 <= Z; Z0++)
            for (len_t ir = 0; ir < Nr; ir++, ionOffset++)
                ni[idx+Z0*Nr+ir] = dynamic_densities[ionOffset+ir];
    }

    eqsys->SetInitialValue(OptionConstants::UQTY_ION_SPECIES, ni, t0);
}

