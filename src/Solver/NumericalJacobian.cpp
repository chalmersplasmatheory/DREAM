/**
 * Routines for numerically evaluating the jacobian numerically
 * in the non-linear solver.
 */

#include <iostream>
#include "DREAM/Solver/SolverNonLinear.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/UnknownQuantity.hpp"


using namespace DREAM;


/**
 * Evaluate jacobian numerically using finite difference
 * of the residual function.
 *
 * jac: Jacobian matrix.
 */
void SolverNonLinear::_EvaluateJacobianNumerically(
    FVM::BlockMatrix *jac
) {
    printf("Evaluating Jacobian numerically...   0.00%%");

    len_t nSize = this->unknowns->GetLongVectorSize(this->nontrivial_unknowns);
    const real_t *iniVec = this->unknowns->GetLongVector(this->nontrivial_unknowns);

    real_t *dFVec   = new real_t[nSize];
    real_t *FVec    = new real_t[nSize];
    real_t *iniFVec = new real_t[nSize];
    real_t *xhVec   = new real_t[nSize];
    
    // Copy initial vector to shifted solution vector
    for (len_t i = 0; i < nSize; i++)
        xhVec[i] = iniVec[i];

    jac->SetOffset(0, 0);
    jac->Zero();

    this->_EvaluateF(iniVec, iniFVec, jac);

    const real_t h = 1e-6, hDefault = 10;
    len_t index = 0;
    for (auto uid : this->nontrivial_unknowns) {
        FVM::UnknownQuantity *uqn = this->unknowns->GetUnknown(uid);
        const len_t N = uqn->NumberOfElements();

        // Differentiate w.r.t. index...
        for (len_t _i = 0; _i < N; _i++, index++) {
            // Restore previous element
            if (index > 0)
                xhVec[index-1] = iniVec[index-1];

            // Determine derivative step length
            real_t hStep;
            if (iniVec[index] == 0)
                hStep = hDefault;
            else
                hStep = h*iniVec[index];

            xhVec[index] += hStep;

            // Evaluate F(x+h)
            this->_EvaluateF(xhVec, FVec, jac);

            // Evaluate dF/dx_index
            for (len_t j = 0; j < nSize; j++)
                dFVec[j] = (FVec[j]-iniFVec[j]) / hStep;

            // Set Jacobian column
            for (len_t j = 0; j < nSize; j++) {
                if (dFVec[j] != 0)
                    jac->SetElement(j, index, dFVec[j], INSERT_VALUES);
            }

            printf("\b\b\b\b\b\b\b%6.2f%%", double(index)/double(nSize-1)*100);
            std::cout << std::flush;
        }
    }

    printf("\n");

    jac->Assemble();

    delete [] xhVec;
    delete [] iniFVec;
    delete [] FVec;
    delete [] dFVec;
}

/**
 * Evaluate the non-linear function 'F'.
 *
 * xVec: Point in which to evaluate the function.
 * FVec: Contains function value on return.
 * jac:  Associated jacobian (for obtaining vector/matrix structure).
 */
void SolverNonLinear::_EvaluateF(const real_t *xVec, real_t *FVec, FVM::BlockMatrix *jac) {
    len_t offset = 0;
    for (auto uqnid : this->nontrivial_unknowns) {
        unknowns->Store(uqnid, xVec, offset);
        offset += unknowns->GetUnknown(uqnid)->NumberOfElements();
    }

    this->RebuildTerms(this->CurrentTime(), this->CurrentTimeStep());
    this->BuildVector(this->CurrentTime(), this->CurrentTimeStep(), FVec, jac);
}

