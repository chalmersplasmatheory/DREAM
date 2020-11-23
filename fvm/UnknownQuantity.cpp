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
 *
 * sf:       SFile object to use for writing the file.
 * path:     Path in the SFile to save data to.
 * saveMeta: If true, also saves grid data for the quantity.
 * current:  If true, saves only data for the most recent iteration/time step.
 */
void UnknownQuantity::SaveSFile(SFile *sf, const string& path, bool saveMeta) {
    this->data->SaveSFile(sf, this->name, path, "", saveMeta);

    string var = path + "/" + this->name;

    sf->WriteAttribute_string(var, "description", this->description);
    sf->WriteAttribute_string(var, "equation", this->description_eqn);
}

/**
 * Save this unknown quantity to the given SFile.
 *
 * sf:       SFile object to use for writing the file.
 * path:     Path in the SFile to save data to.
 * saveMeta: If true, also saves grid data for the quantity.
 * current:  If true, saves only data for the most recent iteration/time step.
 */
void UnknownQuantity::SaveSFileCurrent(SFile *sf, const string& path, bool saveMeta) {
    this->data->SaveSFileCurrent(sf, this->name, path, "", saveMeta);

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

