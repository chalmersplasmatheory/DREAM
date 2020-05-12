
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
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
 * Evaluate this equation.
 *
 * vec:      Vector to store evaluated data in.
 * unknowns: List of unknowns.
 */
void UnknownQuantityEquation::Evaluate(real_t *vec, FVM::UnknownQuantityHandler *unknowns) {
    for (auto it = equations.begin(); it != equations.end(); it++) {
        FVM::UnknownQuantity *uqty = unknowns->GetUnknown(it->first);
        it->second->Evaluate(vec, uqty->GetData());
    }
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
 * If this quantity is predetermined, this routine
 * returns the predetermined parameter. Otherwise,
 * this routine returns a nullptr.
 */
FVM::PredeterminedParameter *UnknownQuantityEquation::GetPredetermined() {
    if (!IsPredetermined())
        return nullptr;
    else
        return equations.begin()->second->GetPredetermined();
}

/**
 * Returns 'true' if this unknown quantity is
 * "evaluable", i.e. can be evaluated without solving
 * an equation.
 */
bool UnknownQuantityEquation::IsEvaluable() {
    bool ev = true;
    for (auto it = equations.begin(); it != equations.end(); it++)
        ev = (ev &&  it->second->IsEvaluable());

    return ev;
}

/**
 * Returns 'true' if this unknown quantity is
 * "predetermined", i.e. set by some external time evolution.
 */
bool UnknownQuantityEquation::IsPredetermined() {
    // A predetermined quantity has an equation that is
    // 
    //   x = x_0(t)
    //
    // where 'x_0(t)' denotes the predetermined value at time t.
    // Hence, if there are multiple equations (= applied to
    // different unknowns) this is not a predetermined quantity.
    if (equations.size() == 1)
        return (equations.begin()->second->IsPredetermined());
    else
        return false;
}

void UnknownQuantityEquation::RebuildEquations(
    const real_t t, const real_t dt, FVM::UnknownQuantityHandler *uqty
) {
    for (auto it = equations.begin(); it != equations.end(); it++)
        it->second->RebuildTerms(t, dt, uqty);
}

/**
 * Evaluate the function vector of this equation.
 *
 * vec:      Vector to store evaluated function in.
 * unknowns: List of unknowns 
 * jac:      Associated jacobian matrix.
 */
void UnknownQuantityEquation::SetVectorElements(
    real_t *vec, FVM::UnknownQuantityHandler *unknowns,
    FVM::BlockMatrix *jac
) {
    for (auto it = equations.begin(); it != equations.end(); it++) {
        FVM::UnknownQuantity *uqty = unknowns->GetUnknown(it->first);
        PetscInt vecoffs = jac->GetOffset(it->first);
        it->second->SetVectorElements(vec + vecoffs, uqty->GetData());
    }
}
