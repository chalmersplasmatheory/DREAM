/**
 * Implementation of a general diffusion term.
 */

#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;

/**
 * Constructor.
 *
 * rg:                Grid on which to define this diffusion term.
 * allocCoefficients: If 'true', allocates memory for the diffusion
 *                    coefficients owned by this term. If 'false',
 *                    it is expected that the memory to use for the
 *                    diffusion coefficients is set using a call to
 *                    'SetCoefficients()' immediately after creation.
 */
DiffusionTerm::DiffusionTerm(Grid *rg, bool allocCoefficients)
    : EquationTerm(rg) {

    if (allocCoefficients)
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

    // XXX Here we explicitly assume that n1[i] = n1[i+1]
    // at all radii
    this->drr[nr] = new real_t[n1[nr-1]*n2[nr-1]];

    this->coefficientsShared = false;
}

/**
 * Deallocates the memory used by the diffusion coefficients.
 */
void DiffusionTerm::DeallocateCoefficients() {
    if (drr != nullptr) {
        for (len_t i = 0; i < grid->GetNr()+1; i++)
            delete [] drr[i];
        delete [] drr;
    }
    if (d11 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d11[i];
        delete [] d11;
    }
    if (d12 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d12[i];
        delete [] d12;
    }
    if (d21 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d21[i];
        delete [] d21;
    }
    if (d22 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] d22[i];
        delete [] d22;
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
    DeallocateCoefficients();

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
 * rhs: Right-hand-side of equation (not side).
 */
void DiffusionTerm::SetMatrixElements(Matrix *mat, real_t*) {
    #define f(K,I,J,V) mat->SetElement(offset+j*np1+i, offset + ((K)-ir)*np2*np1 + (J)*np1 + (I), (V))
    #   include "DiffusionTerm.set.cpp"
    #undef f
}
/*void DiffusionTerm::SetMatrixElements(Matrix *mat, real_t*) {
    const len_t nr = grid->GetNr();
    len_t offset = 0;

    const real_t
        *dr   = grid->GetRadialGrid()->GetDr(),
        *dr_f = grid->GetRadialGrid()->GetDr_f();
    
    // Iterate over interior radial grid points
    for (len_t ir = 0; ir < nr; ir++) {
        const MomentumGrid *mg = grid->GetMomentumGrid(ir);

        const len_t
            np1 = mg->GetNp1(),
            np2 = mg->GetNp2();

        const real_t
            *Vp     = grid->GetVp(ir),
            *Vp_fr  = grid->GetVp_fr(ir),
            *Vp_fr1 = grid->GetVp_fr(ir+1),
            *Vp_f1  = grid->GetVp_f1(ir),
            *Vp_f2  = grid->GetVp_f2(ir),
            *dp1    = mg->GetDp1(),
            *dp2    = mg->GetDp2(),
            *dp1_f  = mg->GetDp1_f(),
            *dp2_f  = mg->GetDp2_f();

        for (len_t j = 0; j < np2; j++) {
            for (len_t i = 0; i < np1; i++) {
                real_t S;

                /////////////////////////
                // RADIUS
                /////////////////////////
                
                #define f(K,V) mat->SetElement(offset + j*np1 + i, \
                    offset+((K)-ir)*np1*np2 + j*np1 + i, \
                    (V))
                
                
                // Phi^(r)_{k-1/2}
                if (ir > 0) {
                    // XXX: Here, we explicitly assume that the momentum grids are
                    // the same at all radii, so that p at (ir, i, j) = p at (ir+1, i, j)
                    S = Drr(ir, i, j)*Vp_fr[j*np1+i] / (dr[ir]*dr_f[ir-1]*Vp[j*np1+i]);
                    f(ir-1, +S);
                    f(ir,   -S);
                }

                // Phi^(r)_{k+1/2}
                if (ir < nr-1) {
                    // XXX: Here, we explicitly assume that the momentum grids are
                    // the same at all radii, so that p at (ir, i, j) = p at (ir+1, i, j)
                    S = Drr(ir+1, i, j)*Vp_fr1[j*np1+i] / (dr[ir]*dr_f[ir]*Vp[j*np1+i]);
                    f(ir,   -S);
                    f(ir+1, +S);
                }

                #undef f
                
                #define f(I,J,V) mat->SetElement(offset+j*np1+i, offset + ((J)*np1) + (I), (V))
                /////////////////////////
                // MOMENTUM 1/1
                /////////////////////////
                // Phi^(1)_{i-1/2,j}
                if (i > 0) {
                    S = D11(ir, i, j)*Vp_f1[j*(np1+1)+i] / (dp1[i]*dp1_f[i-1]*Vp[j*np1+i]);
                    f(i-1, j, +S);
                    f(i, j,   -S);
                }

                // Phi^(1)_{i+1/2,j}
                if (i < np1-1) {
                    S = D11(ir, i+1, j)*Vp_f1[j*(np1+1)+i+1] / (dp1[i]*dp1_f[i]*Vp[j*np1+i]);
                    f(i+1, j, +S);
                    f(i,   j, -S);
                }
                
                /////////////////////////
                // MOMENTUM 2/2
                /////////////////////////
                // Phi^(2)_{i-1/2,j}
                if (j > 0) {
                    S = D22(ir, i, j)*Vp_f2[j*np1+i] / (dp2[j]*dp2_f[j-1]*Vp[j*np1+i]);
                    f(i, j,   -S);
                    f(i, j-1, +S);
                }

                // Phi^(2)_{i+1/2,j}
                if (j < np2-1) {
                    S = D22(ir, i, j)*Vp_f2[(j+1)*np1+i] / (dp2[j]*dp2_f[j]*Vp[j*np1+i]);
                    f(i, j+1, +S);
                    f(i, j,   -S);
                }
                
                /////////////////////////
                // MOMENTUM 1/2
                /////////////////////////
                // Phi^(1)_{i-1/2,j}
                if (i > 0 && (j > 0 && j < np2-1)) {
                    S = D12(ir, i, j)*Vp_f1[j*(np1+1)+i] / (dp1[i]*(dp2_f[j]+dp2_f[j-1])*Vp[j*np1+i]);
                    f(i,   j+1, -S);
                    f(i-1, j+1, -S);
                    f(i,   j-1, +S);
                    f(i-1, j-1, +S);
                }

                // Phi^(1)_{i+1/2,j}
                if (i < np1-1 && (j > 0 && j < np2-1)) {
                    S = D12(ir,i+1,j)*Vp_f1[j*(np1+1)+i+1]/(dp1[i]*(dp2_f[j]+dp2_f[j-1])*Vp[j*np1+i]);
                    f(i+1, j+1, +S);
                    f(i,   j+1, +S);
                    f(i+1, j-1, -S);
                    f(i,   j-1, -S);
                }
                
                /////////////////////////
                // MOMENTUM 2/1
                /////////////////////////
                // Phi^(2)_{i,j-1/2}
                if (j > 0 && (i > 0 && i < np1-1)) {
                    S = D21(ir,i,j)*Vp_f2[j*np1+i] / (dp2[j]*(dp1_f[i]+dp1_f[i-1])*Vp[j*np1+i]);
                    f(i+1, j-1, -S);
                    f(i+1, j,   -S);
                    f(i-1, j-1, +S);
                    f(i-1, j,   +S);
                }

                // Phi^(2)_{i,j+1/2}
                if (j < np2-1 && (i > 0 && i < np1-1)) {
                    S = D21(ir,i,j+1)*Vp_f2[(j+1)*np1+i] / (dp2[j]*(dp1_f[i]+dp1_f[i-1])*Vp[j*np1+i]);
                    f(i+1, j+1, +S);
                    f(i+1, j,   +S);
                    f(i-1, j+1, -S);
                    f(i-1, j,   -S);
                }
                
                
                #undef f
            }
        }

        offset += np1*np2;
    }
}*/

/**
 * Instead of building a linear operator (matrix) to apply to a vector
 * 'x', this routine builds immediately the resulting vector.
 *
 * vec: Vector to set elements of.
 * x:   Input x vector.
 */
void DiffusionTerm::SetVectorElements(real_t *vec, const real_t *x) {
    #define f(K,I,J,V) vec[offset+j*np1+i] += (V)*x[offset+((K)-ir)*np2*np1 + (J)*np1 + (I)]
    #   include "DiffusionTerm.set.cpp"
    #undef f
}

