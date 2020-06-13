/**
 * Implementation of the NIST ADS data interface.
 */

#include <unordered_map>
#include "DREAM/DREAMException.hpp"
#include "DREAM/NIST.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
NIST::NIST() {
    // Add total binding energies
    for (len_t i = 0; i < nist_binding_n; i++) {
        struct nist_data *b = nist_binding_table+i;
        this->binding[b->Z] = b;
    }

    // Add ionization energies
    for (len_t i = 0; i < nist_ionization_n; i++) {
        struct nist_data *b = nist_ionization_table+i;
        this->ionization[b->Z] = b;
    }
}


/**
 * Retrieve the total binding energy for the ion with
 * atomic charge Z, in charge state Z0.
 *
 * Z:  Atomic charge of ion.
 * Z0: Charge state of ion.
 */
real_t NIST::GetBindingEnergy(const len_t Z, const len_t Z0) const {
    return this->_GetEnergy(Z, Z0, this->binding);
}

/**
 * Retrieve the ionization energy for the ion with
 * atomic charge Z, in charge state Z0.
 *
 * Z:  Atomic charge of ion.
 * Z0: Charge state of ion.
 */
real_t NIST::GetIonizationEnergy(const len_t Z, const len_t Z0) const {
    return this->_GetEnergy(Z, Z0, this->ionization);
}

/**
 * General, internal method for retrieving NIST data.
 */
real_t NIST::_GetEnergy(
    const len_t Z, const len_t Z0,
    const unordered_map<len_t, struct nist_data*>& arr
) const {
    auto it = arr.find(Z);
    if (it == arr.end())
        throw DREAMException(
            "No binding energy available for ions with charge Z = "
            LEN_T_PRINTF_FMT ".",
            Z
        );

    if (Z0 == Z)
        return 0;
    else
        return it->second->data[Z0];
}

