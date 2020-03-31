/**
 * Implementation of a general advection term.
 */

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM::FVM;

/**
 * Constructor.
 */
AdvectionTerm::AdvectionTerm(Grid *rg)
    : EquationTerm(rg) {
    
    this->AllocateCoefficients();
}

/**
 * Destructor.
 */
AdvectionTerm::~AdvectionTerm() {
    if (!this->coefficientsShared)
        DeallocateCoefficients();
}

/**
 * Allocate new memory for the advection coefficients
 * based on the current grid sizes.
 */
void AdvectionTerm::AllocateCoefficients() {
    if (!this->coefficientsShared)
        DeallocateCoefficients();

    this->fr = new real_t*[nr+1];
    this->f1 = new real_t*[nr];
    this->f2 = new real_t*[nr];

    for (len_t i = 0; i < nr; i++) {
        this->fr[i] = new real_t[n1[i]*n2[i]];
        this->f1[i] = new real_t[(n1[i]+1)*n2[i]];
        this->f2[i] = new real_t[n1[i]*(n2[i]+1)];
    }

    // TODO What about this point???
    //this->fr[nr] = new real_t[???];
    // XXX: Here we assume that the momentum grid is the same
    // at all radial grid points, so that n1_{nr+1/2} = n1_{nr-1/2}
    // (and the same for n2)
    this->fr[nr] = new real_t[n1[nr-1]*n2[nr-1]];

    this->coefficientsShared = false;
}

/**
 * Deallocates the memory used by the advection coefficients.
 */
void AdvectionTerm::DeallocateCoefficients() {
    if (f2 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] f2[i];

        delete [] f2;
    }
    if (f1 != nullptr) {
        for (len_t i = 0; i < grid->GetNr(); i++)
            delete [] f1[i];

        delete [] f1;
    }
    if (fr != nullptr) {
        for (len_t i = 0; i < grid->GetNr()+1; i++)
            delete [] fr[i];

        delete [] fr;
    }
}

/**
 * Assign the memory regions to store the coefficients
 * of this term. This means that we will assume that the
 * memory region is 'shared', and will leave it to someone
 * else to 'free' the memory later on.
 *
 * fr: List of radial advection coefficients.
 * f1: List of first momentum advection coefficients.
 * f2: List of second momentum advection coefficients.
 */
void AdvectionTerm::SetCoefficients(real_t **fr, real_t **f1, real_t **f2) {
    this->fr = fr;
    this->f1 = f1;
    this->f2 = f2;

    this->coefficientsShared = true;
}

/**
 * This function is called whenever the computational grid is
 * re-built, in case the grid has been re-sized (in which case we might
 * need to re-allocate memory for the advection coefficients)
 */
bool AdvectionTerm::GridRebuilt() {
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
void AdvectionTerm::SetMatrixElements(Matrix *mat) {

    const len_t nr = grid->GetNr();
    len_t offset = 0;

    const real_t
        *dr = grid->GetRadialGrid()->GetDr();

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
            *dp2    = mg->GetDp2();

        /*const real_t
            *Fr  = fr[ir],
            *Fr1 = fr[ir+1],
            *F1  = f1[ir],
            *F2  = f2[ir];*/

        for (len_t j = 0; j < np2; j++) {
            // Evaluate flux in first point
            //real_t S = F1[j*(np1+1) + 1] * h2_f1[j*(np1+1) + 1] * h3_f1[j*(np1+1) + 1] / dp1[1];

            for (len_t i = 0; i < np1; i++) {
                /////////////////////////
                // RADIUS
                /////////////////////////
                #define f(K,V) mat->SetElement(offset + j*np1 + i, offset + ((K)-ir)*np2*np1 + j*np1 + i, (V))
                // XXX: Here we assume that the momentum grid is the same at all
                // radial points.
                //
                // OTHERWISE...:
                // This term is pretty difficult, since we need to evaluate the
                // coefficients in different radial grid points, which in general
                // have different momentum grids and thus require interpolation
                // across grids. For this application, we should be able to assume
                // that momentum grids at all radii use the same coordinate systems.
                // In general, it is much more difficult, though (so it would require
                // a bit more thinking if we wanted to interpolate generally between
                // two different momentum grids)

                // Phi^(r)_{ir-1/2,i,j}
                if (ir > 0) {
                    real_t S = Fr(ir, i, j) * Vp_fr[j*np1+i] / (Vp[j*np1+i] * dr[ir]);
                    f(ir-1, -S * (1-deltar[ir][j*np1 + i]));
                    f(ir,   -S * deltar[ir][j*np1 + i]);
                }

                // Phi^(r)_{ir+1/2,i,j}
                if (ir < nr-1) {
                    real_t S = Fr(ir+1, i, j) * Vp_fr1[j*np1+i] / (Vp[j*np1+i] * dr[ir]);
                    f(ir,   S * (1-deltar[ir+1][j*np1 + i]));
                    f(ir+1, S * deltar[ir+1][j*np1 + i]);
                }
                
                #undef f
                
                /////////////////////////
                // MOMENTUM 1
                /////////////////////////
                #define f(I,J,V) mat->SetElement(offset + j*np1 + i, offset + ((J)*np1) + (I), (V))

                // Phi^(1)_{i-1/2,j}
                if (i > 0) {
                    real_t S = F1(ir, i, j) * Vp_f1[j*(np1+1) + i] / (Vp[j*np1+i]*dp1[i]);
                    f(i-1, j,-S * (1-delta1[ir][j*np1 + i]));
                    f(i,   j,-S * delta1[ir][j*np1 + i]);
                }

                // Phi^(1)_{i+1/2,j}
                if (i < np1-1) {
                    real_t S = F1(ir, i+1, j) * Vp_f1[j*(np1+1) + i+1] / (Vp[j*np1+i]*dp1[i]);
                    f(i,   j, S * (1-delta1[ir][j*np1 + i+1]));
                    f(i+1, j, S * delta1[ir][j*np1 + i+1]);
                }

                /////////////////////////
                // MOMENTUM 2
                /////////////////////////
                // Phi^(2)_{i,j-1/2}
                if (j > 0) {
                    real_t S2m = F2(ir, i, j) * Vp_f2[j*np1+i] / (Vp[j*np1+i]*dp2[j]);
                    f(i, j-1,-S2m * (1-delta2[ir][j*np1+i]));
                    f(i, j,  -S2m * delta2[ir][j*np1+i]);
                }

                // Phi^(2)_{i,j+1/2}
                if (j < np2-1) {
                    real_t S2p = F2(ir, i, j+1) * Vp_f2[(j+1)*np1+i] / (Vp[j*np1+i]*dp2[j+1]);
                    f(i, j,   S2p * (1-delta2[ir][(j+1)*np1+i]));
                    f(i, j+1, S2p * delta2[ir][(j+1)*np1+i]);
                }

                #undef f
            }
        }

        offset += np1*np2;
    }
}

