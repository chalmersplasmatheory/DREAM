/**
 * Implementation of ion equations.
 */

#include "DREAM/Settings/SimulationGenerator.hpp"


using namespace DREAM;
using namespace std;


#define MODULENAME "ions"

/**
 * Define options for the ions.
 *
 * set: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_Ions(Settings *s) {
    const len_t dims[1] = {0};

    s->DefineSetting(MODULENAME "/Z", "List of atomic charge numbers", 1, dims, (int_t*)nullptr);
    s->DefineSetting(MODULENAME "/types", "Method to use for determining ion charge distributions", 1, dims, (int_t*)nullptr);

    //DefineDataIonR(MODULENAME, s, "density");
    DefineDataIonRT(MODULENAME, s, "prescribed");
}

/**
 * Construct the equation governing the evolution of the
 * ion densities.
 */
void SimulationGenerator::ConstructEquation_Ions(EquationSystem *eqsys, Settings *s) {
    Grid *fluidGrid = eqsys->GetFluidGrid();

    len_t nZ, ntypes;
    int_t *Z  = s->GetIntegerArray(MODULENAME "/Z", 1, &nZ);
    int_t *itypes = s->GetIntegerArray(MODULENAME "/types", 1, &ntypes);

    // Verify that exactly one type per ion species is given
    if (nZ != ntypes)
        throw SettingsException(
            "ions: Expected the lengths of 'Z' and 'types' to match."
        );

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
                prescribed_indices[nZ_prescribed++]++;
                break;

            case OptionConstants::ION_DATA_TYPE_DYNAMIC:
            case OptionConstants::ION_DATA_TYPE_EQUILIBRIUM:
                dynamic_indices[nZ_dynamic++]++;
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
        MODULENAME, fluidGrid->GetRadialGrid(), s, nZ_dynamic, "densities"
    );
    IonInterpolator1D *prescribed_densities = LoadDataIonRT(
        MODULENAME, fluidGrid->GetRadialGrid(), s, nZ0_prescribed, "prescribed"
    );

    IonHandler *ih = new IonHandler(fluidGrid->GetRadialGrid(), unknowns, Z, nZ);
    eqsys->SetIonHandler(ih);

    // Initialize ion equations
    Equation *eqn = new Equation(fluidGrid);

    IonPrescribedParameter *ipp = nullptr;
    if (nZ0_prescribed > 0)
        ipp = new IonPrescribedParameter(fluidGrid, nZ_prescribed, prescribed_indices, prescribed_densities);

    // Construct dynamic equations
    for (len_t iZ = 0; iZ < nZ; iZ++) {
        switch (types[iZ]) {
            case OptionConstants::ION_DATA_PRESCRIBED: break;

            case OptionConstants::ION_DATA_TYPE_DYNAMIC:
                eqn->AddTerm(new IonRateEquation(
                    fluidGrid, ih, iZ, eqsys->GetADAS(), unknowns
                ));
                break;

            case OptionConstants::ION_DATA_TYPE_EQUILIBRIUM:
                throw NotImplementedException(
                    "The 'equilibrium' ion equation type has not been implemented yet."
                );

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
}

