/**
 * Implementation of a general diffusion term.
 */

#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/RadialGrid.hpp"


using namespace TQS::FVM;

/**
 * Constructor.
 */
DiffusionTerm::DiffusionTerm(RadialGrid *rg)
    : EquationTerm(rg) {

    this->AllocateCoefficients();
}

/**
 * Destructor.
 */
DiffusionTerm::~DiffusionTerm() {
    if (!this->coefficientsShared)
        DeallocateCoefficients();
}

/**
 * Allocate new memory for the diffusion coefficients.
 */
void DiffusionTerm::AllocateCoefficients() {
    if (!this->coefficientsShared)
        DeallocateCoefficients();
    
    this->drr = new real_t*[nr+1];
    this->d11 = new real_t*[nr];
    this->d12 = new real_t*[nr];
    this->d21 = new real_t*[nr];
    this->d22 = new real_t*[nr];

    for (len_t i = 0; i < nr; i++) {
        this->drr[i] = new real_t[n1[i]*n2[i]];
        this->d11[i] = new real_t[(n1[i]+1)*n2[i]];
        this->d12[i] = new real_t[(n1[i]+1)*n2[i]];
        this->d21[i] = new real_t[n1[i]*(n2[i]+1)];
        this->d22[i] = new real_t[n1[i]*(n2[i]+1)];
    }

    this->coefficientsShared = false;
}

/**
 * Deallocates the memory used by the diffusion coefficients.
 */
void DiffusionTerm::DeallocateCoefficients() {
    if (drr != nullptr) {
        for (len_t i = 0; i < grid->GetNr()+1; i++)
            delete [] drr[i];
    }
    if (d11 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d11[i];
    }
    if (d12 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d12[i];
    }
    if (d21 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d21[i];
    }
    if (d22 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d22[i];
    }
}

/**
 * Assign the memory regions to store the coefficients
 * of this term. This means that we will assume that the
 * memory region is 'shared', and will leave it to someone
 * else to 'free' the memory later on.
 *
 * drr: List of R/R diffusion coefficients.
 * d11: List of p1/p1 diffusion coefficients.
 * d12: List of p1/p2 diffusion coefficients.
 * d21: List of p2/p1 diffusion coefficients.
 * d22: List of p2/p2 diffusion coefficients.
 */
void DiffusionTerm::SetCoefficients(
    real_t **drr,
    real_t **d11, real_t **d12,
    real_t **d21, real_t **d22
) {
    this->drr = drr;
    this->d11 = d11;
    this->d12 = d12;
    this->d21 = d21;
    this->d22 = d22;

    this->coefficientsShared = true;
}

/**
 * This function is called whenever the computational grid is
 * re-built, in case the grid has been re-sized (in which case
 * we might need to re-allocate memory for the diffusion coefficients)
 */
bool DiffusionTerm::GridRebuilt() {
    this->EquationTerm::GridRebuilt();

    // Do not re-build if our coefficients are owned by someone else
    if (this->coefficientsShared)
        return false;

    this->AllocateCoefficients();

    return true;
}

/**
 * Build the matrix elements for this operator.
 *
 * mat: Matrix to build elements of.
 */
void DiffusionTerm::SetMatrixElements(Matrix *mat) {
    const len_t nr = grid->GetNr();
    len_t offset = 0;
    
    // Iterate over interior radial grid points
    for (len_t ir = 0; ir < nr; ir++) {
        const MomentumGrid *mg = grid->GetMomentumGrid(ir);

        const len_t
            np1 = mg->GetNp1(),
            np2 = mg->GetNp2();

        const real_t
            *h2_f1 = mg->GetH2_f1(),
            *h3_f1 = mg->GetH3_f1(),
            *h1_f2 = mg->GetH1_f2(),
            *h3_f2 = mg->GetH3_f2(),
            *dp1   = mg->GetDp1(),
            *dp2   = mg->GetDp2(),
            *dp1_f = mg->GetDp1_f(),
            *dp2_f = mg->GetDp2_f();

        const real_t
            *D11 = d11[ir],
            *D12 = d12[ir],
            *D21 = d21[ir],
            *D22 = d22[ir];

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                /////////////////////////
                // RADIUS
                /////////////////////////
                
                #define f(I,J,V) mat->SetElement(offset+j*np1+i, offset + ((J)*np1) + (I), (V))
                /////////////////////////
                // MOMENTUM 1/1
                /////////////////////////
                
                /////////////////////////
                // MOMENTUM 1/2
                /////////////////////////
                
                /////////////////////////
                // MOMENTUM 2/1
                /////////////////////////
                
                /////////////////////////
                // MOMENTUM 2/2
                /////////////////////////
                
                #undef f
            }
        }

        offset += np1*np2;
    }
}

