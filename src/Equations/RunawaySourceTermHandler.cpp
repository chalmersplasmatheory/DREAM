/**
 * A collection of fluid runaway source terms.
 */

#include <algorithm>
#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/RunawaySourceTermHandler.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
RunawaySourceTermHandler::RunawaySourceTermHandler() { }

/**
 * Destructor.
 */
RunawaySourceTermHandler::~RunawaySourceTermHandler() {
    this->applyToAll([this](FVM::EquationTerm *eqterm) {
        delete eqterm;
    });
}

/**
 * Apply the given operator to all source terms that
 * have been initialized.
 */
void RunawaySourceTermHandler::applyToAll(std::function<void(FVM::EquationTerm*)> op) {
    if (this->avalanche != nullptr)
        op(this->avalanche);
    if (this->compton != nullptr)
        op(this->compton);
    if (this->dreicer != nullptr)
        op(this->dreicer);
}
void RunawaySourceTermHandler::applyToAll(const std::function<void(FVM::EquationTerm*)> op) const {
    if (this->avalanche != nullptr)
        op(this->avalanche);
    if (this->compton != nullptr)
        op(this->compton);
    if (this->dreicer != nullptr)
        op(this->dreicer);
}

/**
 * Method called when the simulation grid is rebuilt.
 */
/*bool RunawaySourceTermHandler::GridRebuilt() {
    bool rebuilt = this->EquationTerm::GridRebuilt();

    // Rebuild all sources
    this->applyToAll([this,&rebuilt](FVM::EquationTerm *src) {
        rebuilt = rebuilt || src->GridRebuilt();
    });

    return rebuilt;
}*/

/**
 * Get number of non-zero elements per row for the
 * matrix and jacobian.
 */
/*len_t RunawaySourceTermHandler::GetNumberOfNonZerosPerRow() const {
    len_t nnz = 0;

    this->applyToAll([this,&nnz](FVM::EquationTerm *eq) {
        nnz = max(nnz, eq->GetNumberOfNonZerosPerRow());
    });

    return nnz;
}*/

/**
 * Get number of non-zero elements per row for the
 * matrix and jacobian.
 */
/*len_t RunawaySourceTermHandler::GetNumberOfNonZerosPerRow_jac() const {
    len_t nnz = 0;

    this->applyToAll([this,&nnz](FVM::EquationTerm *eq) {
        nnz = max(nnz, eq->GetNumberOfNonZerosPerRow_jac());
    });

    return nnz;
}*/

/**
 * Rebuild all source tems.
 */
/*void RunawaySourceTermHandler::Rebuild(
    const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns
) {
    this->applyToAll([this,&t,&dt,&unknowns](FVM::EquationTerm *eq) {
        eq->Rebuild(t, dt, unknowns);
    });
}*/

/**
 * Set the Jacobian elements for all source terms.
 */
/*void RunawaySourceTermHandler::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *x
) {
    this->applyToAll([this,&uqtyId,&derivId,&jac,&x](FVM::EquationTerm *eq) {
        eq->SetJacobianBlock(uqtyId, derivId, jac, x);
    });
}*/

/**
 * Set the linear operator matrix elements for all source terms.
 */
/*void RunawaySourceTermHandler::SetMatrixElements(
    FVM::Matrix *mat, real_t *rhs
) {
    this->applyToAll([this,&mat,&rhs](FVM::EquationTerm *eq) {
        eq->SetMatrixElements(mat, rhs);
    });
}*/

/**
 * Set the Jacobian elements for all source terms.
 */
/*void RunawaySourceTermHandler::SetVectorElements(
    real_t *vec, const real_t *x
) {
    this->applyToAll([this,&vec,&x](FVM::EquationTerm *eq) {
        eq->SetVectorElements(vec, x);
    });
}*/

/**
 * Add the available source terms to the given operators.
 *
 * op_nRE:  Operator operating on n_re.
 * op_nTot: Operator operating on n_tot.
 */
void RunawaySourceTermHandler::AddToOperators(
    FVM::Operator *op_nRE, FVM::Operator *op_nTot
) {
    // n_re
    if (this->avalanche != nullptr) {
        if (op_nRE == nullptr)
            throw DREAMException(
                "RunawaySourceTermHandler: Avalanche enabled, but no operator "
                "for n_re provided."
            );
        else
            op_nRE->AddTerm(this->avalanche);
    }
    if (this->dreicer != nullptr) {
        if (op_nRE == nullptr)
            throw DREAMException(
                "RunawaySourceTermHandler: Avalanche enabled, but no operator "
                "for n_re provided."
            );
        else
            op_nRE->AddTerm(this->dreicer);
    }

    // n_tot
    if (this->compton != nullptr) {
        if (op_nTot == nullptr)
            throw DREAMException(
                "RunawaySourceTermHandler: Avalanche enabled, but no operator "
                "for n_tot provided."
            );
        else
            op_nTot->AddTerm(this->compton);
    }
}

