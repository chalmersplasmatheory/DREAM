
#include <map>
#include "DREAM/DREAMException.hpp"
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
 * uqtyId:   ID of unknown quantity to evaluate.
 * vec:      Vector to store evaluated data in.
 * unknowns: List of unknowns.
 */
void UnknownQuantityEquation::Evaluate(
    const len_t uqtyId, real_t *vec, FVM::UnknownQuantityHandler *unknowns
) {
    FVM::Operator *eqn = nullptr;

    for (auto it = equations.begin(); it != equations.end(); it++) {
        FVM::UnknownQuantity *uqty = unknowns->GetUnknown(it->first);
        
        // This equation will be used later to "solve" for
        // the unknown quantity, so store it for now...
        if (it->first == uqtyId) {
            eqn = it->second;

            if (!eqn->IsEvaluable())
                throw DREAMException(
                    "%s: quantity equation is not evaluable because the diagonal term is not evaluable.",
                    uqty->GetName().c_str()
                );
        } else
            it->second->Evaluate(vec, uqty->GetData());
    }

    // Predetermined equation terms have an implicit identity term
    // applied to them.
    if (!IsPredetermined()) {
        if (eqn == nullptr)
            throw DREAMException(
                "%s: quantity equation is not evaluable because it has no diagonal term.",
                unknowns->GetUnknown(uqtyId)->GetName().c_str()
            );

        // Apply transformation which solves for the unknown quantity
        eqn->EvaluableTransform(vec);
    } else
        eqn->Evaluate(vec, nullptr);
}

/**
 * Returns true if this equation contains a transient term.
 */
bool UnknownQuantityEquation::HasTransientTerm() const {
	for (auto it = equations.begin(); it != equations.end(); it++)
		if (it->second->HasTransientTerm())
			return true;
	
	return false;
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
    /*bool ev = true;
    for (auto it = equations.begin(); it != equations.end(); it++)
        ev = (ev &&  it->second->IsEvaluable());*/

    if (equations.find(this->uqtyId) == equations.end())
        return false;

    FVM::Operator *eqn = equations[this->uqtyId];
    return eqn->IsEvaluable();
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
 * uqn_id:   ID of unknown quantity.
 * vec:      Vector to store evaluated function in.
 * unknowns: List of unknowns 
 * jac:      Associated jacobian matrix.
 */
void UnknownQuantityEquation::SetVectorElements(
    real_t *vec, FVM::UnknownQuantityHandler *unknowns
) {
    for (auto it = equations.begin(); it != equations.end(); it++) {
        FVM::UnknownQuantity *uqty = unknowns->GetUnknown(it->first);
        it->second->SetVectorElements(vec, uqty->GetData());
    }
}
