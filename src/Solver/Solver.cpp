/**
 * Implementation of common routines for the 'Solver' routines.
 */

#include <vector>
#include "DREAM/Solver/Solver.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/UnknownQuantity.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
Solver::Solver(
    FVM::UnknownQuantityHandler *unknowns,
    vector<UnknownQuantityEquation*> *unknown_equations
)
    : unknowns(unknowns), unknown_equations(unknown_equations) {
}


/**
 * Build a jacobian matrix for the equation system.
 *
 * t:    Time to build the jacobian matrix for.
 * dt:   Length of time step to take.
 * mat:  Matrix to use for storing the jacobian.
 */
void Solver::BuildJacobian(const real_t, const real_t, FVM::BlockMatrix *jac) {
    // Reset jacobian matrix
    jac->Zero();

    // Iterate over (non-trivial) unknowns (i.e. those which appear
    // in the matrix system)
    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        len_t uqn_id = nontrivial_unknowns[i];
        UnknownQuantityEquation *eqn = unknown_equations->at(uqn_id);
        const real_t *x = unknowns->GetUnknownData(uqn_id);
        
        // Iterate over each equation
        for (auto it = eqn->GetEquations().begin(); it != eqn->GetEquations().end(); it++) {
            // "Differentiate with respect to the unknowns which
            // appear in the matrix"
            //   d (eqn_it) / d x_j
            for (len_t j = 0; j < nontrivial_unknowns.size(); j++) {
                len_t derivId = nontrivial_unknowns[j];
                jac->SelectSubEquation(this->unknownToMatrixMapping[derivId], this->unknownToMatrixMapping[it->first]);
                it->second->SetJacobianBlock(derivId, it->first, jac, x);
            }
        }
    }

    jac->PartialAssemble();

    // Apply boundary conditions which overwrite elements
    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        len_t uqn_id = nontrivial_unknowns[i];
        UnknownQuantityEquation *eqn = unknown_equations->at(uqn_id);
        const real_t *x = unknowns->GetUnknownData(uqn_id);
        
        // Iterate over each equation
        for (auto it = eqn->GetEquations().begin(); it != eqn->GetEquations().end(); it++) {
            // "Differentiate with respect to the unknowns which
            // appear in the matrix"
            //   d (eqn_it) / d x_j
            for (len_t j = 0; j < nontrivial_unknowns.size(); j++) {
                len_t derivId = nontrivial_unknowns[j];
                jac->SelectSubEquation(this->unknownToMatrixMapping[derivId], this->unknownToMatrixMapping[it->first]);
                it->second->SetJacobianBlockBC(derivId, it->first, jac, x);
            }
        }
    }

    jac->Assemble();
}

/**
 * Build a linear operator matrix for the equation system.
 *
 * t:    Time to build the jacobian matrix for.
 * dt:   Length of time step to take.
 * mat:  Matrix to use for storing the jacobian.
 * rhs:  Right-hand-side in equation.
 */
void Solver::BuildMatrix(const real_t, const real_t, FVM::BlockMatrix *mat, real_t *S) {
    // Reset matrix and rhs
    mat->Zero();
    for (len_t i = 0; i < matrix_size; i++)
        S[i] = 0;

    // Build matrix
    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        len_t uqn_id = nontrivial_unknowns[i];
        UnknownQuantityEquation *eqn = unknown_equations->at(uqn_id);

        for (auto it = eqn->GetEquations().begin(); it != eqn->GetEquations().end(); it++) {
            mat->SelectSubEquation(this->unknownToMatrixMapping[uqn_id], this->unknownToMatrixMapping[it->first]);
            PetscInt vecoffs = mat->GetOffset(this->unknownToMatrixMapping[it->first]);
            it->second->SetMatrixElements(mat, S + vecoffs);
        }
    }

    mat->Assemble();
}

/**
 * Build a function vector for the equation system.
 *
 * t:   Time to build the function vector for.
 * dt:  Length of time step to take.
 * vec: Vector to store evaluated equations in.
 * jac: Associated jacobian matrix.
 */
void Solver::BuildVector(const real_t, const real_t, real_t *vec, FVM::BlockMatrix *jac) {
    // Reset function vector
    for (len_t i = 0; i < matrix_size; i++)
        vec[i] = 0;

    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        unknown_equations->at(i)->SetVectorElements(vec, unknowns, jac);
    }
}

/**
 * Initialize this solver.
 *
 * size:     Number of elements in full unknown vector.
 *           (==> jacobian is of size 'size-by-size').
 * unknowns: List of indices of unknowns to include in the
 *           function vectors/matrices.
 */
void Solver::Initialize(const len_t size, vector<len_t>& unknowns) {
    this->matrix_size = size;

    // Copy list of non-trivial unknowns (those which will
    // appear in the matrices that are built later on)
    nontrivial_unknowns = unknowns;

    this->initialize_internal(size, unknowns);
}

/**
 * Rebuild all equation terms in the equation system for
 * the specified time.
 *
 * t:  Time for which to rebuild the equation system.
 * dt: Length of time step to take next.
 */
void Solver::RebuildTerms(const real_t t, const real_t dt) {
    // Update prescribed quantities
    const len_t N = unknowns->Size();
    for (len_t i = 0; i < N; i++) {
        FVM::UnknownQuantity *uqty = unknowns->GetUnknown(i);
        UnknownQuantityEquation *eqn = unknown_equations->at(i);

        if (eqn->IsPredetermined()) {
            eqn->RebuildEquations(t, dt, unknowns);
            FVM::PredeterminedParameter *pp = eqn->GetPredetermined();
            uqty->Store(pp->GetData(), 0, true);
        }
    }

    for (len_t i = 0; i < nontrivial_unknowns.size(); i++) {
        len_t uqnId = nontrivial_unknowns[i];
        UnknownQuantityEquation *eqn = unknown_equations->at(uqnId);

        for (auto it = eqn->GetEquations().begin(); it != eqn->GetEquations().end(); it++) {
            it->second->RebuildTerms(t, dt, unknowns);
        }
    }
}

