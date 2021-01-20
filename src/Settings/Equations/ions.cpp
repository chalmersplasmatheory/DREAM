/**
 * Implementation of ion equations.
 */

#include "DREAM/Equations/Fluid/IonPrescribedParameter.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include "DREAM/Equations/Fluid/IonKineticIonizationTerm.hpp"
#include "DREAM/Equations/Fluid/IonTransientTerm.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/Operator.hpp"


using namespace DREAM;
using namespace std;


#define MODULENAME "eqsys/n_i"

/**
 * Define options for the ions.
 *
 * s: Settings object to define options in.
 */
void SimulationGenerator::DefineOptions_Ions(Settings *s) {
    const len_t dims[1] = {0};

    s->DefineSetting(MODULENAME "/names", "Names of each ion species", (const string)"");
    s->DefineSetting(MODULENAME "/Z", "List of atomic charge numbers", 1, dims, (int_t*)nullptr);
    s->DefineSetting(MODULENAME "/types", "Method to use for determining ion charge distributions", 1, dims, (int_t*)nullptr);
    s->DefineSetting(MODULENAME "/tritiumnames", "Names of the tritium ion species", (const string)"");
    s->DefineSetting(MODULENAME "/ionization", "Model to use for ionization", (int_t) OptionConstants::EQTERM_IONIZATION_MODE_FLUID);
    s->DefineSetting(MODULENAME "/typeTi", "Model to use for ion heat equation", (int_t) OptionConstants::UQTY_T_I_NEGLECT);

    DefineDataIonR(MODULENAME, s, "initial");
    DefineDataIonR(MODULENAME, s, "initialTi");
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
 * Returns the number of ion species in the system
 * (i.e. the number of elements "divided by nr (number of radial points)"
 * in the "NI_DENS" and "WI_ENER" unknown quantities)
 *
 * set: Settings object to load charge state number from.
 */
len_t SimulationGenerator::GetNumberOfIonSpecies(Settings *s) {
    len_t nZ;
    s->GetIntegerArray(MODULENAME "/Z", 1, &nZ, false);
    return nZ;
}

/**
 * Construct the equation governing the evolution of the
 * ion densities for each charge state.
 */
void SimulationGenerator::ConstructEquation_Ions(EquationSystem *eqsys, Settings *s, ADAS *adas) {
    const real_t t0 = 0;
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    len_t nZ, ntypes;
    const int_t *_Z  = s->GetIntegerArray(MODULENAME "/Z", 1, &nZ);
    const int_t *itypes = s->GetIntegerArray(MODULENAME "/types", 1, &ntypes);

    // Parse list of ion names (stored as one contiguous string,
    // each substring separated by ';')
    vector<string> ionNames = s->GetStringList(MODULENAME "/names");

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

    // Get list of tritium species
    vector<string> tritiumNames = s->GetStringList(MODULENAME "/tritiumnames");

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

    // Verify that all non-prescribed elements are in ADAS
    for (len_t i = 0; i < nZ; i++) {
        if (!adas->HasElement(Z[i]) && types[i] != OptionConstants::ION_DATA_PRESCRIBED)
            throw SettingsException(
                "ions: The DREAM ADAS database does not contain '%s' (Z = " LEN_T_PRINTF_FMT ")",
                ionNames[i].c_str(), Z[i]
            );
    }

    /////////////////////
    /// LOAD ION DATA ///
    /////////////////////
    // Count number of prescribed/dynamic charge states
    len_t nZ0_prescribed=0, nZ_prescribed=0, nZ_dynamic=0, nZ0_dynamic=0;
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
                nZ0_dynamic += Z[i] + 1;
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
        MODULENAME, fluidGrid->GetRadialGrid(), s, nZ0_dynamic, "initial"
    );
    MultiInterpolator1D *prescribed_densities = LoadDataIonRT(
        MODULENAME, fluidGrid->GetRadialGrid(), s, nZ0_prescribed, "prescribed"
    );

    IonHandler *ih = new IonHandler(fluidGrid->GetRadialGrid(), eqsys->GetUnknownHandler(), Z, nZ, ionNames, tritiumNames);
    eqsys->SetIonHandler(ih);

    // Initialize ion equations
    FVM::Operator *eqn = new FVM::Operator(fluidGrid);

    OptionConstants::eqterm_ionization_mode ionization_mode = 
        (enum OptionConstants::eqterm_ionization_mode)s->GetInteger(MODULENAME "/ionization");
    FVM::Operator *Op_kiniz = nullptr; 
    FVM::Operator *Op_kiniz_re = nullptr; 
    // TODO: simplify the bool logic below if possible
    bool includeKineticIonization = (ionization_mode == OptionConstants::EQTERM_IONIZATION_MODE_KINETIC) 
                                 || (ionization_mode == OptionConstants::EQTERM_IONIZATION_MODE_KINETIC_APPROX_JAC);
    if(includeKineticIonization && !(eqsys->HasHotTailGrid()||eqsys->HasRunawayGrid()))
        throw SettingsException("Invalid ionization mode: cannot use kinetic ionization without a kinetic grid.");
    bool collfreqModeIsFull = (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode")
        == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL;
    bool addFluidIonization = !(includeKineticIonization && eqsys->HasHotTailGrid() && collfreqModeIsFull);
    bool addFluidJacobian = (includeKineticIonization && eqsys->HasHotTailGrid() && (ionization_mode==OptionConstants::EQTERM_IONIZATION_MODE_KINETIC_APPROX_JAC));
    IonPrescribedParameter *ipp = nullptr;
    if (nZ0_prescribed > 0)
        ipp = new IonPrescribedParameter(fluidGrid, ih, nZ_prescribed, prescribed_indices, prescribed_densities);

    const len_t id_ni = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    // Construct dynamic equations
    len_t nDynamic = 0, nEquil = 0;
    for (len_t iZ = 0; iZ < nZ; iZ++) {
        switch (types[iZ]) {
            case OptionConstants::ION_DATA_PRESCRIBED: 
                break;
            // 'Dynamic' and 'Equilibrium' differ by a transient term
            case OptionConstants::ION_DATA_TYPE_DYNAMIC:
                nDynamic++;
                eqn->AddTerm(
                    new IonTransientTerm(fluidGrid, ih, iZ, id_ni)
                );
                [[fallthrough]];
            case OptionConstants::ION_DATA_EQUILIBRIUM:
                nEquil++;
                eqn->AddTerm(new IonRateEquation(
                    fluidGrid, ih, iZ, adas, eqsys->GetUnknownHandler(),
                    addFluidIonization, addFluidJacobian
                ));
                if(includeKineticIonization){
                    if(eqsys->HasHotTailGrid()) { // add kinetic ionization to hot-tail grid
                        if(Op_kiniz == nullptr)
                            Op_kiniz = new FVM::Operator(eqsys->GetHotTailGrid());
                        Op_kiniz->AddTerm(new IonKineticIonizationTerm(
                            fluidGrid, eqsys->GetHotTailGrid(), id_ni, 
                            eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT), eqsys->GetUnknownHandler(), 
                            ih, iZ, ionization_mode, eqsys->GetHotTailGridType()==OptionConstants::MOMENTUMGRID_TYPE_PXI, 
                            collfreqModeIsFull, eqsys->GetUnknownID(OptionConstants::UQTY_F_HOT)
                        ));
                    }
                    // TODO: always include RE ionization (as long as HasRunawayGrid), but
                    //       consider using a simple jacobian (assume Ion_re ~ n_re)
                    if(eqsys->HasRunawayGrid()) {
                        if(Op_kiniz_re == nullptr)
                            Op_kiniz_re = new FVM::Operator(eqsys->GetRunawayGrid());
                        Op_kiniz_re->AddTerm(new IonKineticIonizationTerm(
                            fluidGrid, eqsys->GetRunawayGrid(), id_ni, 
                            eqsys->GetUnknownID(OptionConstants::UQTY_F_RE), eqsys->GetUnknownHandler(), 
                            ih, iZ, ionization_mode, eqsys->GetRunawayGridType()==OptionConstants::MOMENTUMGRID_TYPE_PXI, 
                            false, eqsys->GetUnknownID(OptionConstants::UQTY_F_RE)
                        )); 
                    }
                }
                break;

            default:
                throw SettingsException(
                    "ions: Unrecognized ion model type specified: %d.",
                    types[iZ]
                );
        }
    }

    // Set equation description
    string desc;
    if (ipp != nullptr && nEquil > 0) {
        if (nEquil == nDynamic)
            desc = "Prescribed + dynamic";
        else
            desc = "Prescribed + dynamic + equilibrium";
    } else if (ipp != nullptr)
        desc = "Fully prescribed";
    else {
        if (nEquil == nDynamic)
            desc = "Fully dynamic";
        else if (nDynamic == 0)
            desc = "Fully equilibrium";
        else
            desc = "Dynamic + equilibrium";
    }
    if((nEquil>0) && includeKineticIonization)
        desc += " (kinetic ionization)";

    if (ipp != nullptr)
        eqn->AddTerm(ipp);

    eqsys->SetOperator(id_ni, id_ni, eqn, desc);
    if(Op_kiniz != nullptr)
        eqsys->SetOperator(id_ni, OptionConstants::UQTY_F_HOT, Op_kiniz, desc);
    if(Op_kiniz_re != nullptr)
        eqsys->SetOperator(id_ni, OptionConstants::UQTY_F_RE, Op_kiniz_re, desc);

    // Initialize dynamic ions
    const len_t Nr = fluidGrid->GetNr();
    real_t *ni = new real_t[ih->GetNzs() * Nr];

    for (len_t i = 0; i < ih->GetNzs() * Nr; i++)
        ni[i] = 0;

    // Begin by evaluating prescribed densities
    if (ipp != nullptr) {
        ipp->Rebuild(t0, 1, nullptr);
        ipp->Evaluate(ni);
    }

    // ...and then fill in with the initial dynamic ion values
    for (len_t i = 0, ionOffset = 0; i < nZ_dynamic; i++) {
        len_t Z   = ih->GetZ(dynamic_indices[i]);
        len_t idx = ih->GetIndex(dynamic_indices[i], 0);

        for (len_t Z0 = 0; Z0 <= Z; Z0++) {
            for (len_t ir = 0; ir < Nr; ir++)
                ni[(idx+Z0)*Nr+ir] = dynamic_densities[ionOffset+ir];
            ionOffset += Nr;
        }
    }

    eqsys->SetInitialValue(id_ni, ni, t0);
    ih->Rebuild();

    delete [] types;
}
