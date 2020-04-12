/**
 * Implementation of the FVM UnknownQuantity class.
 */

#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/UnknownQuantity.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

using namespace DREAM::FVM;
using namespace std;


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

