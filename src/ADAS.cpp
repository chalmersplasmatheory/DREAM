/**
 * Implementation of the object which handles ADAS rate coefficients and
 * generates the necessary interpolation objects.
 */

#include <map>
#include "DREAM/ADAS.hpp"
#include "DREAM/adasdata.h"
#include "DREAM/IO.hpp"


using namespace std;
using namespace DREAM;


const len_t
    ADAS::IDX_ACD=0,
    ADAS::IDX_CCD=1,
    ADAS::IDX_SCD=2,
    ADAS::IDX_PLT=3,
    ADAS::IDX_PRB=4;


/**
 * Constructor.
 */
ADAS::ADAS(const gsl_interp2d_type *interp) {
    // Generate interpolator objects
    for (len_t i = 0; i < adas_rate_n; i++) {
        struct adas_rate *ar = (adas_rate_table+i);

        ADASRateInterpolator **ari = new ADASRateInterpolator*[4];

        #define INITADAS(type,shiftZ0) \
            new ADASRateInterpolator( \
                ar->Z, ar-> type ## _nn, ar-> type ## _nT, \
                ar-> type ## _n, ar-> type ## _T, \
                ar-> type , shiftZ0, interp \
            )

        ari[IDX_ACD] = INITADAS(acd, true);
        ari[IDX_CCD] = INITADAS(ccd, false);
        ari[IDX_SCD] = INITADAS(scd, false);
        ari[IDX_PLT] = INITADAS(plt, false);
        ari[IDX_PRB] = INITADAS(prb, true);

        intp[ar->Z] = ari;

        #undef INITADAS
    }
}


/**
 * Destructor.
 */
ADAS::~ADAS() {
    for (auto it = intp.begin(); it != intp.end(); it++) {
        // Iterate over data types
        for (len_t i = 0; i < 4; i++)
            delete it->second[i];

        delete [] it->second;
    }
}


/**
 * Check if the element with the specified atomic charge
 * exists in this ADAS database.
 */
bool ADAS::HasElement(const len_t Z) const {
    return (intp.find(Z) != intp.end());
}


/**
 * (private)
 * Returns an iterator to the entry in the interpolator list
 * for the element with the specified charge.
 *
 * Z: Charge of element to get iterator to.
 */
map<len_t, ADASRateInterpolator**>::const_iterator ADAS::get_element(const len_t Z) const {
    map<len_t, ADASRateInterpolator**>::const_iterator it = intp.find(Z);
    if (it == intp.end())
        throw ADASException(
            "Element with charge '" LEN_T_PRINTF_FMT "' not in DREAM ADAS database.",
            Z
        );

    return it;
}

/**
 * Getters for ADAS data.
 */
ADASRateInterpolator *ADAS::GetACD(const len_t Z) const {
    return get_element(Z)->second[IDX_ACD];
}
ADASRateInterpolator *ADAS::GetCCD(const len_t Z) const {
    return get_element(Z)->second[IDX_CCD];
}
ADASRateInterpolator *ADAS::GetSCD(const len_t Z) const {
    return get_element(Z)->second[IDX_SCD];
}
ADASRateInterpolator *ADAS::GetPLT(const len_t Z) const {
    return get_element(Z)->second[IDX_PLT];
}
ADASRateInterpolator *ADAS::GetPRB(const len_t Z) const {
    return get_element(Z)->second[IDX_PRB];
}

/**
 * Print a list of all elements available in the DREAM ADAS database.
 */
void ADAS::PrintElements() const {
    DREAM::IO::PrintInfo("ELEMENTS IN ADAS");
    for (len_t i = 0; i < adas_rate_n; i++) {
        struct adas_rate *ar = (adas_rate_table+i);

        DREAM::IO::PrintInfo("  (%2" LEN_T_PRINTF_FMT_PART ") %s", ar->Z, ar->name);
    }
}

