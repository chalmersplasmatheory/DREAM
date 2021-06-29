/**
 * Implementation of the object which handles ADAS rate coefficients and
 * generates the necessary interpolation objects.
 */

#include <unordered_map>
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


#define N_ADAS_RATES 5

/**
 * Constructor.
 */
ADAS::ADAS(const gsl_interp2d_type *interp) {
    // Generate interpolator objects
    for (len_t i = 0; i < adas_rate_n; i++) {
        struct adas_rate *ar = (adas_rate_table+i);

        ADASRateInterpolator **ari = new ADASRateInterpolator*[N_ADAS_RATES];

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

        len_t idx = get_isotope_index(ar->Z, ar->A);
        intp[idx] = ari;

        #undef INITADAS
    }
}


/**
 * Destructor.
 */
ADAS::~ADAS() {
    for (auto it = intp.begin(); it != intp.end(); it++) {
        // Iterate over data types
        for (len_t i = 0; i < N_ADAS_RATES; i++)
            delete it->second[i];

        delete [] it->second;
    }
}


/**
 * (private)
 * Returns the index into the ADAS rate array for the specified
 * isotope.
 */
len_t ADAS::get_isotope_index(const len_t Z, const len_t A) const {
    return (Z*MAX_ATOMIC_MASS + (A==0 ? Z : A));
}


/**
 * Check if the element with the specified atomic charge
 * exists in this ADAS database.
 */
bool ADAS::HasElement(const len_t Z, const len_t A) const {
    len_t idx = get_isotope_index(Z, A);
    return (intp.find(idx) != intp.end());
}


/**
 * (private)
 * Returns an iterator to the entry in the interpolator list
 * for the element with the specified charge.
 *
 * Z: Charge of element to get iterator to.
 */
unordered_map<len_t, ADASRateInterpolator**>::const_iterator ADAS::get_element(const len_t Z, const len_t A) const {
    len_t idx = get_isotope_index(Z, A);

    unordered_map<len_t, ADASRateInterpolator**>::const_iterator it = intp.find(idx);
    if (it == intp.end()) {
        if (A == 0)
            throw ADASException(
                "Element with charge Z=" LEN_T_PRINTF_FMT " not in DREAM ADAS database.",
                Z
            );
        else
            throw ADASException(
                "Element with charge Z=" LEN_T_PRINTF_FMT " and mass A=" LEN_T_PRINTF_FMT " not in DREAM ADAS database.",
                Z, A
            );
    }

    return it;
}

/**
 * Getters for ADAS data.
 */
ADASRateInterpolator *ADAS::GetACD(const len_t Z, const len_t A) const {
    return get_element(Z, A)->second[IDX_ACD];
}
ADASRateInterpolator *ADAS::GetCCD(const len_t Z, const len_t A) const {
    return get_element(Z, A)->second[IDX_CCD];
}
ADASRateInterpolator *ADAS::GetSCD(const len_t Z, const len_t A) const {
    return get_element(Z, A)->second[IDX_SCD];
}
ADASRateInterpolator *ADAS::GetPLT(const len_t Z, const len_t A) const {
    return get_element(Z, A)->second[IDX_PLT];
}
ADASRateInterpolator *ADAS::GetPRB(const len_t Z, const len_t A) const {
    return get_element(Z, A)->second[IDX_PRB];
}

/**
 * Print a list of all elements available in the DREAM ADAS database.
 */
void ADAS::PrintElements() const {
    DREAM::IO::PrintInfo("ELEMENTS IN ADAS");
    for (len_t i = 0; i < adas_rate_n; i++) {
        struct adas_rate *ar = (adas_rate_table+i);

        DREAM::IO::PrintInfo(
            "  (A=%2" LEN_T_PRINTF_FMT_PART ", Z=%2" LEN_T_PRINTF_FMT_PART ") %s",
            ar->A, ar->Z, ar->name
        );
    }
}

