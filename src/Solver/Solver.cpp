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
    // in the matrix system), corresponding to blocks in F and
    // rows in the Jacobian matrix.
    for (len_t uqnId : nontrivial_unknowns) {
        UnknownQuantityEquation *eqn = unknown_equations->at(uqnId);
        const real_t *x = unknowns->GetUnknownData(uqnId);
        map<len_t, len_t>& utmm = this->unknownToMatrixMapping;
        
        // Iterate over each equation term
        for (auto it = eqn->GetEquations().begin(); it != eqn->GetEquations().end(); it++) {
            // If the unknown quantity to which this operator is applied is
            // trivial (and thus not part of the matrix system), it's derivative
            // does not appear in the Jacobian (and is most likely 0 anyway), and
            // so we silently skip it
            if (utmm.find(it->first) == utmm.end())
                continue;

            // "Differentiate with respect to the unknowns which
            // appear in the matrix"
            //   d (F_uqnId) / d x_derivId
            for (len_t derivId : nontrivial_unknowns) {
                jac->SelectSubEquation(utmm[uqnId], utmm[derivId]);

                // - in the equation for                           x_uqnId
                // - differentiate the operator that is applied to x_it
                // - with respect to                               x_derivId
                it->second->SetJacobianBlock(it->first, derivId, jac, x);
            }
        }
    }

    jac->PartialAssemble();

    // Apply boundary conditions which overwrite elements
    for (len_t uqnId : nontrivial_unknowns) {
        UnknownQuantityEquation *eqn = unknown_equations->at(uqnId);
        const real_t *x = unknowns->GetUnknownData(uqnId);
        map<len_t, len_t>& utmm = this->unknownToMatrixMapping;
        
        // Iterate over each equation
        for (auto it = eqn->GetEquations().begin(); it != eqn->GetEquations().end(); it++) {
            // Skip trivial unknowns
            if (utmm.find(it->first) == utmm.end())
                continue;

            // "Differentiate with respect to the unknowns which
            // appear in the matrix"
            //   d (eqn_uqnId) / d x_derivId
            for (len_t derivId : nontrivial_unknowns) {
                jac->SelectSubEquation(this->unknownToMatrixMapping[uqnId], this->unknownToMatrixMapping[derivId]);
                // For logic, see comment in the for-loop above
                it->second->SetJacobianBlockBC(it->first, derivId, jac, x);
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
    for (len_t uqnId : nontrivial_unknowns) {
        UnknownQuantityEquation *eqn = unknown_equations->at(uqnId);
        map<len_t, len_t>& utmm = this->unknownToMatrixMapping;

        for (auto it = eqn->GetEquations().begin(); it != eqn->GetEquations().end(); it++) {
            if (utmm.find(it->first) != utmm.end()) {
                mat->SelectSubEquation(utmm[uqnId], utmm[it->first]);
                PetscInt vecoffs = mat->GetOffset(utmm[uqnId]);
                it->second->SetMatrixElements(mat, S + vecoffs);

            // The unknown to which this operator should be applied is a
            // "trivial" unknown quantity, meaning it does not appear in the
            // equation system matrix. We therefore build it as part of the
            // RHS vector.
            } else {
                PetscInt vecoffs = mat->GetOffset(utmm[uqnId]);
                const real_t *data = unknowns->GetUnknownData(it->first);
                it->second->SetVectorElements(S + vecoffs, data);
            }
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
        len_t uqn_id = nontrivial_unknowns[i];
        len_t vecoffset = jac->GetOffset(i);
        unknown_equations->at(uqn_id)->SetVectorElements(vec+vecoffset, unknowns);
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
    // Rebuild collision handlers and RunawayFluid
    if (this->cqh_hottail != nullptr)
        this->cqh_hottail->Rebuild();
    if (this->cqh_runaway != nullptr)
        this->cqh_runaway->Rebuild();


    bool useApproximateEceffMethod = true; 
    /**
     * true:  Use approximate pitch distribution which has 
     *        an error <10% in the exponent. ~3 times less 
     *        CPU use than 'false' setting.
     * false: Evaluate pitch distribution via gsl integration 
     */
    this->REFluid -> Rebuild(useApproximateEceffMethod);

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

