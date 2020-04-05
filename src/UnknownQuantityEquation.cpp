
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

using namespace DREAM;


/**
 * Destructor.
 */
UnknownQuantityEquation::~UnknownQuantityEquation() {
    for (auto it = equations.begin(); it != equations.end(); it++)
        delete it->second;
}

/**
 * Returns the number of non-zero elements per row
 * in the linear operator matrix constructed for this
 * unknown.
 * (NOTE: this is just an estimate)
 */
len_t UnknownQuantityEquation::NumberOfNonZeros() {
    len_t nnz = 0;
    for (auto it = equations.begin(); it != equations.end(); it++)
        nnz += it->second->GetNumberOfNonZerosPerRow();

    return nnz;
}

/**
 * Returns the number of non-zero elements per row
 * in the jacobian matrix constructed for this unknown.
 * (NOTE: this is just an estimate)
 */
len_t UnknownQuantityEquation::NumberOfNonZeros_jac() {
    len_t nnz = 0;
    for (auto it = equations.begin(); it != equations.end(); it++)
        nnz += it->second->GetNumberOfNonZerosPerRow_jac();

    return nnz;
}

/**
 * If this quantity is prescribed, this routine
 * returns the prescribed parameter. Otherwise,
 * this routine returns a nullptr.
 */
FVM::PrescribedParameter *UnknownQuantityEquation::GetPrescribed() {
    if (!IsPrescribed())
        return nullptr;
    else
        return equations.begin()->second->GetPrescribed();
}


/**
 * Returns 'true' if this unknown quantity is
 * "prescribed", i.e. set by some external time evolution.
 */
bool UnknownQuantityEquation::IsPrescribed() {
    // A prescribed quantity has an equation that is
    // 
    //   x = x_0(t)
    //
    // where 'x_0(t)' denotes the prescribed value at time t.
    // Hence, if there are multiple equations (= applied to
    // different unknowns) this is not prescribed quantity.
    if (equations.size() == 1)
        return (equations.begin()->second->IsPrescribed());
    else
        return false;
}

void UnknownQuantityEquation::RebuildEquations(
    const real_t t, const real_t dt, FVM::UnknownQuantityHandler *uqty
) {
    for (auto it = equations.begin(); it != equations.end(); it++)
        it->second->RebuildTerms(t, dt, uqty);
}

