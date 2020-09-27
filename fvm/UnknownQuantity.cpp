/**
 * Implementation of the FVM UnknownQuantity class.
 */

#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/UnknownQuantity.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

using namespace DREAM::FVM;
using namespace std;


/**
 * Save this unknown quantity to the given SFile.
 */
void UnknownQuantity::SaveSFile(SFile *sf, const string& path, bool saveMeta) {
    this->data->SaveSFile(sf, this->name, path, "", saveMeta);

    string var = path + "/" + this->name;

    sf->WriteAttribute_string(var, "description", this->description);
    sf->WriteAttribute_string(var, "equation", this->description_eqn);
}

/**
 * Set the initial value of the specified unknown quantity. If
 * the initial value has previously been specified, it is overwritten.
 *
 * val: Initial value of the quantity.
 * t0:  Initial time.
 */
void UnknownQuantity::SetInitialValue(const real_t *val, const real_t t0) {
    this->data->SetInitialValue(val, t0);
}

