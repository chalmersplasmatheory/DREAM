
#include <map>
#include "DREAM/DREAMException.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

using namespace DREAM;


/**
 * Constructor.
 */
UnknownQuantityEquation::UnknownQuantityEquation(
	len_t uqtyId, FVM::UnknownQuantity *uqty, const std::string& desc,
	EquationTriggerCondition *condition
) :uqtyId(uqtyId), uqty(uqty), condition(condition), description(desc) {
	
	if (condition == nullptr) {
		const len_t N = uqty->NumberOfElements();
		this->eqn_cache = new real_t[N];
		this->nElements = N;
	}
}

/**
 * Destructor.
 */
UnknownQuantityEquation::~UnknownQuantityEquation() {
    for (auto it = equations.begin(); it != equations.end(); it++)
        delete it->second;
	
	for (auto it = equations_alt.begin(); it != equations_alt.end(); it++)
		delete it->second;
	
	if (this->eqn_cache != nullptr)
		delete [] this->eqn_cache;
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
    FVM::Operator *eqn = nullptr, *eqn_a = nullptr;

    for (auto it = this->equations.begin(); it != this->equations.end(); it++) {
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

	// Evaluate alternative equation
	if (this->condition != nullptr) {
		for (auto it = this->equations.begin(); it != this->equations.end(); it++) {
			FVM::UnknownQuantity *uqty = unknowns->GetUnknown(it->first);
			
			// This equation will be used later to "solve" for
			// the unknown quantity, so store it for now...
			if (it->first == uqtyId) {
				eqn_a = it->second;

				if (!eqn_a->IsEvaluable())
					throw DREAMException(
						"%s: quantity equation is not evaluable because the diagonal term is not evaluable.",
						uqty->GetName().c_str()
					);
			} else
				it->second->Evaluate(eqn_cache, uqty->GetData());
		}
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
    } else {
        eqn->Evaluate(vec, nullptr);

		if (eqn_a != nullptr) {
			eqn_a->Evaluate(eqn_cache, nullptr);

			// Set the elements for which the condition is triggered
			for (len_t i = 0; i < nElements; i++) {
				if (this->condition->IsTriggered(i))
					vec[i] = eqn_cache[i];
			}
		}
	}
}


/**
 * Returns true if this equation contains a transient term.
 */
bool UnknownQuantityEquation::HasTransientTerm() const {
	for (auto it = this->equations.begin(); it != this->equations.end(); it++)
		if (it->second->HasTransientTerm())
			return true;
	
	if (this->condition != nullptr)
		for (auto it = this->equations_alt.begin(); it != this->equations_alt.end(); it++)
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
    for (auto it = this->equations.begin(); it != this->equations.end(); it++)
        nnz += it->second->GetNumberOfNonZerosPerRow();
	
	if (this->condition != nullptr) {
		len_t nnz_a = 0;
		for (auto it = this->equations_alt.begin(); it != this->equations_alt.end(); it++)
			nnz_a += it->second->GetNumberOfNonZerosPerRow();

		if (nnz_a > nnz)
			nnz = nnz_a;
	}

    return nnz;
}

/**
 * Returns the number of non-zero elements per row
 * in the jacobian matrix constructed for this unknown.
 * (NOTE: this is just an estimate)
 */
len_t UnknownQuantityEquation::NumberOfNonZeros_jac() {
    len_t nnz = 0;
    for (auto it = this->equations.begin(); it != this->equations.end(); it++)
        nnz += it->second->GetNumberOfNonZerosPerRow_jac();

	if (this->condition != nullptr) {
		len_t nnz_a = 0;
		for (auto it = this->equations_alt.begin(); it != this->equations_alt.end(); it++)
			nnz_a += it->second->GetNumberOfNonZerosPerRow_jac();

		if (nnz_a > nnz)
			nnz = nnz_a;
	}

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
	// Check that there is at least one term applied to 'uqtyId'
	// (an evaluable term must contain an identity term, applied
	//  applied to 'uqtyId')
    if (this->equations.find(this->uqtyId) == this->equations.end())
        return false;

    FVM::Operator *eqn = this->equations[this->uqtyId];
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
    if (equations.size() == 1) {
		bool eq0 = equations.begin()->second->IsPredetermined();
		if (equations_alt.size() == 0)
			return eq0;
		else
			return false;
    } else
        return false;
}

void UnknownQuantityEquation::RebuildEquations(
    const real_t t, const real_t dt, FVM::UnknownQuantityHandler *uqty
) {
	if (this->condition != nullptr)
		this->condition->CheckCondition(t, uqty);

	for (auto it = this->equations.begin(); it != this->equations.end(); it++)
		it->second->RebuildTerms(t, dt, uqty);
	
	// If available, also rebuild alternative equation
	if (this->condition != nullptr)
		for (auto it = this->equations_alt.begin(); it != this->equations_alt.end(); it++)
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
    for (auto it = this->equations.begin(); it != this->equations.end(); it++) {
        FVM::UnknownQuantity *uqty = unknowns->GetUnknown(it->first);
        it->second->SetVectorElements(vec, uqty->GetData());
    }

	if (this->condition != nullptr) {
		for (len_t i = 0; i < nElements; i++)
			this->eqn_cache[i] = 0;

		for (auto it = this->equations_alt.begin(); it != this->equations_alt.end(); it++) {
			FVM::UnknownQuantity *uqty = unknowns->GetUnknown(it->first);
			it->second->SetVectorElements(eqn_cache, uqty->GetData());
		}

		// Set the elements for which the condition is triggered
		const len_t N = this->condition->GetNCells();
		for (len_t i = 0; i < N; i++)
			if (this->condition->IsTriggered(i))
				vec[i] = eqn_cache[i];
	}
}


/**
 * Set the trigger condition to use for switching between
 * different equations.
 */
void UnknownQuantityEquation::SetTriggerCondition(
	EquationTriggerCondition *condition
) {
	this->condition = condition;

	if (condition == nullptr) {
		const len_t N = uqty->NumberOfElements();
		this->eqn_cache = new real_t[N];
		this->nElements = N;
	}
}

