/**
 * A test for the IonRateEquation term in DREAM.
 */

#include <vector>
#include <string>
#include "IonRateEquation.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonRateEquation.hpp"
#include "DREAM/Settings/OptionConstants.hpp"


using namespace DREAMTESTS::_DREAM;
using namespace std;


const len_t N_IONS = 4;
const len_t Z_IONS[N_IONS] = {1,4,10,18};
const char ION_NAMES[N_IONS][3] = {"H","Be","Ne","Ar"};


/**
 * Generate a default ion handler.
 */
DREAM::IonHandler *IonRateEquation::GetIonHandler(
    DREAM::FVM::Grid *g, DREAM::FVM::UnknownQuantityHandler *uqh
) {
    vector<string> tritiumNames(0);
    vector<string> names(N_IONS);
    for (len_t i = 0; i < N_IONS; i++)
        names[i] = ION_NAMES[i];

    return new DREAM::IonHandler(
        g->GetRadialGrid(), uqh, Z_IONS, N_IONS, names, tritiumNames
    );
}

/**
 * Generate an unknown quantity handler.
 */
DREAM::FVM::UnknownQuantityHandler *IonRateEquation::GetUnknownHandler(DREAM::FVM::Grid *g) {
    DREAM::FVM::UnknownQuantityHandler *uqh = new DREAM::FVM::UnknownQuantityHandler();

    len_t nZ0 = 0;
    for (len_t i = 0; i < N_IONS; i++)
        nZ0 += Z_IONS[i] + 1;

    this->id_ions = uqh->InsertUnknown(DREAM::OptionConstants::UQTY_ION_SPECIES, "0", g, nZ0);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_COLD, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_HOT, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_N_TOT, "0", g);
    uqh->InsertUnknown(DREAM::OptionConstants::UQTY_T_COLD, "0", g);

    // Set initial values
    const len_t N = nZ0*g->GetNr();
    real_t *nions = new real_t[N];
    len_t ionOffset = 0, rOffset = 0;
    for (len_t iIon = 0; iIon < N_IONS; iIon++) {
        for (len_t Z0 = 0; Z0 <= Z_IONS[iIon]; Z0++, ionOffset++) {
            for (len_t ir = 0; ir < g->GetNr(); ir++, rOffset++)
                nions[rOffset] = (ionOffset+1) * 1e18;
        }
    }
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_ION_SPECIES, nions);

    #define SETVAL(NAME, v) do { \
            for (len_t i = 0; i < g->GetNr(); i++) \
                temp[i] = (v); \
            uqh->SetInitialValue((NAME), temp); \
        } while (false)

    // Set electron quantities
    real_t *temp = new real_t[g->GetNr()];
    SETVAL(DREAM::OptionConstants::UQTY_N_COLD, 1e19);
    SETVAL(DREAM::OptionConstants::UQTY_N_HOT,  5e17);
    SETVAL(DREAM::OptionConstants::UQTY_T_COLD, 1000);

    delete [] nions;
    delete [] temp;
    
    return uqh;
}

/**
 * Check if the ion rate equation conserves particle density.
 */
bool IonRateEquation::CheckConservativity() {
    bool success = true;
    DREAM::FVM::Grid *grid = this->InitializeFluidGrid();
    DREAM::FVM::UnknownQuantityHandler *uqh = GetUnknownHandler(grid);
    DREAM::IonHandler *ih = GetIonHandler(grid, uqh);
    ih->Rebuild();
    DREAM::ADAS *adas = new DREAM::ADAS();
    const len_t Nr = grid->GetNr();

    // Set total electron density
    const real_t *ntot = ih->GetFreePlusBoundElectronDensity();
    uqh->SetInitialValue(DREAM::OptionConstants::UQTY_N_TOT, ntot);

    // Construct solution vector
    len_t nSize = uqh->GetUnknown(this->id_ions)->NumberOfElements();
    real_t *vec = new real_t[nSize];
    real_t *nions = uqh->GetUnknown(this->id_ions)->GetData();
    
    // Clear 'vec'
    for (len_t i = 0; i < nSize; i++)
        vec[i] = 0;

    // Construct equation for each ion species
    DREAM::IonRateEquation *ire[N_IONS];
    for (len_t iIon = 0; iIon < N_IONS; iIon++)
        ire[iIon] = new DREAM::IonRateEquation(grid, ih, iIon, adas, uqh, true, true,false);


    // Check the equation for each ion species
    len_t rOffset = 0;
    for (len_t iIon = 0; iIon < N_IONS; iIon++) {
        ire[iIon]->Rebuild(0, 1, uqh);

        for (len_t Z0 = 0; Z0 <= Z_IONS[iIon]; Z0++, rOffset += Nr)
            ire[iIon]->SetCSVectorElements(vec, nions, iIon, Z0, rOffset);
    }

    // Sum of all elements should vanish
    len_t ionOffset = 0;
    for (len_t iIon = 0; iIon < N_IONS; iIon++) {
        for (len_t ir = 0; ir < Nr; ir++) {
            // Sum over charge states
            real_t s = 0, maxval = 0;
            for (len_t Z0 = 0; Z0 <= Z_IONS[iIon]; Z0++) {
                real_t v = vec[(ionOffset+Z0)*Nr + ir];
                s += v;
                
                if (fabs(v) > maxval)
                    maxval = v;
            }

            if (fabs(s) >= 100 * maxval * std::numeric_limits<real_t>::epsilon()) {
                this->PrintError(
                    "Sum of rate equation over all charge states of ion species "
                    "Z=" LEN_T_PRINTF_FMT " at radius ir=" LEN_T_PRINTF_FMT
                    " does not vanish. Sum = %e",
                    Z_IONS[iIon], ir, s
                );
                success = false;
                goto BREAK_LOOP;
            }
        }

        ionOffset += Z_IONS[iIon]+1;
    }

BREAK_LOOP:

    // Deallocate equation terms
    for (len_t iIon = 0; iIon < N_IONS; iIon++)
        delete ire[iIon];

    delete [] vec;
    delete adas;
    delete ih;
    delete uqh;
    delete grid;

    return success;
}

/**
 * Run this test.
 */
bool IonRateEquation::Run(bool) {
    bool success = true;

    if ((success &= CheckConservativity()))
        this->PrintOK("The ion rate equation conserves density.");
    else
        this->PrintError("The ion rate equation does not conserve density.");

    return success;
}

