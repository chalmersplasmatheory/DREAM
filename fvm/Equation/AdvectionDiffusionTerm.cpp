/**
 * Implementation of a combined advection-diffusion term.
 */

#include <algorithm>
#include "FVM/Equation/AdvectionDiffusionTerm.hpp"


using namespace DREAM::FVM;
using namespace std;

/**
 * Add an advection term to this term.
 */
void AdvectionDiffusionTerm::Add(AdvectionTerm *a) {
    a->SetCoefficients(this->fr, this->f1, this->f2, this->f1pSqAtZero);
    a->SetInterpolationCoefficients(this->deltar, this->delta1, this->delta2, this->deltars, this->delta1s, this->delta2s);
    advectionterms.push_back(a);
}

/**
 * Add a diffusion term to this term.
 */
void AdvectionDiffusionTerm::Add(DiffusionTerm *d) {
    d->SetCoefficients(
        this->drr, this->d11, this->d12, this->d21, this->d22
    );
    diffusionterms.push_back(d);
}

/**
 * Returns the number of non-zero elements per row
 * inserted by this term into a linear operator matrix.
 */
len_t AdvectionDiffusionTerm::GetNumberOfNonZerosPerRow() const {
    len_t nnz = 0;
    for (auto it = advectionterms.begin(); it != advectionterms.end(); it++)
        nnz = max(nnz, (*it)->GetNumberOfNonZerosPerRow());

    for (auto it = diffusionterms.begin(); it != diffusionterms.end(); it++)
        nnz = max(nnz, (*it)->GetNumberOfNonZerosPerRow());

    return nnz;
}

/**
 * Returns the number of non-zero elements per row
 * inserted by this term into a jacobian matrix.
 */
len_t AdvectionDiffusionTerm::GetNumberOfNonZerosPerRow_jac() const {
    len_t nnz = 0;
    for (auto it = advectionterms.begin(); it != advectionterms.end(); it++)
        nnz = max(nnz, (*it)->GetNumberOfNonZerosPerRow_jac());

    for (auto it = diffusionterms.begin(); it != diffusionterms.end(); it++)
        nnz = max(nnz, (*it)->GetNumberOfNonZerosPerRow_jac());

    return nnz;
}

/**
 * Reset the advection and diffusion coefficients of this term.
 */
void AdvectionDiffusionTerm::ResetCoefficients() {
    if (this->advectionterms.size() > 0)
        this->AdvectionTerm::ResetCoefficients();
    if (this->diffusionterms.size() > 0)
        this->DiffusionTerm::ResetCoefficients();

}

/**
 * Rebuild this equation term.
 *
 * t: Simulation time to rebuild term for.
 */
void AdvectionDiffusionTerm::Rebuild(const real_t t, const real_t dt, UnknownQuantityHandler *uqty) {
    this->ResetCoefficients();

    // Rebuild advection-diffusion coefficients
    for (auto it = advectionterms.begin(); it != advectionterms.end(); it++){
        (*it)->Rebuild(t, dt, uqty);
    }

    for (auto it = diffusionterms.begin(); it != diffusionterms.end(); it++){
        (*it)->Rebuild(t, dt, uqty);
    }

    // Rebuild interpolation coefficients
    RebuildInterpolationCoefficients(uqty);
}

/**
 * Rebuild the interpolation coefficients.
 */
void AdvectionDiffusionTerm::RebuildInterpolationCoefficients(UnknownQuantityHandler* /*unknowns*/) {
    ResetInterpolationCoefficients();
    switch (this->interpolationMethod) {
        case AD_INTERP_CENTRED:  
            SetInterpolationCoefficientValues(0.5); 
            break;
        case AD_INTERP_BACKWARD: 
            SetInterpolationCoefficientValues(0); 
            break;
        case AD_INTERP_FORWARD:  
            SetInterpolationCoefficientValues(1); 
            break;
        case AD_INTERP_HYBRID_UPWIND:   
            SetInterpolationCoefficientValuesHybridUpwind(); 
            break;
        case AD_INTERP_QUICK:
            SetInterpolationCoefficientQUICK(deltars, fr, FVM::FLUXGRIDTYPE_RADIAL);
            SetInterpolationCoefficientQUICK(delta1s, f1, FVM::FLUXGRIDTYPE_P1);
            SetInterpolationCoefficientQUICK(delta2s, f2, FVM::FLUXGRIDTYPE_P2);
            break;

        default:
            throw EquationTermException(
                "Unrecognized advection/diffusion interpolation method specified: %d.",
                this->interpolationMethod
            );
    }
    ApplyInterpolationCoefficientBoundaryCondition(deltars, AD_BC_MIRRORED, AD_BC_DIRICHLET, FVM::FLUXGRIDTYPE_RADIAL);
    ApplyInterpolationCoefficientBoundaryCondition(delta1s, AD_BC_MIRRORED, AD_BC_DIRICHLET, FVM::FLUXGRIDTYPE_P1);
    ApplyInterpolationCoefficientBoundaryCondition(delta2s, AD_BC_MIRRORED, AD_BC_MIRRORED, FVM::FLUXGRIDTYPE_P2);
}

void AdvectionDiffusionTerm::ResetInterpolationCoefficients(){

    const len_t nr = this->AdvectionTerm::nr;

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t n2 = this->AdvectionTerm::n2[ir];
        const len_t n1 = this->AdvectionTerm::n1[ir];
        // reset coefficients
        for(len_t i=0; i<(n1+1)*n2; i++)
            for(len_t k=0; k<2*stencil_order; k++)
                delta1s[ir][i][k] = 0;
        for(len_t i=0; i<n1*(n2+1); i++)
            for(len_t k=0; k<2*stencil_order; k++)
                delta2s[ir][i][k] = 0;
    }
    const len_t n2 = this->AdvectionTerm::n2[0]; // XXX: assume same grid at all radii
    const len_t n1 = this->AdvectionTerm::n1[0];
    // reset coefficients
    for(len_t ir=0; ir<nr+1; ir++)
        for(len_t i=0; i<n1*n2; i++)
            for(len_t k=0;k<2*stencil_order;k++)
                deltars[ir][i][k] = 0;
}


/**
 * Sets interpolation coefficients written on "flux limited" form
 * with a (piecewise) linear limiter function psi(r) = a*r + b:
 *      f_{i+1/2} = f_i + 0.5*dx_i*psi(r_{i+1/2})*(f_i - f_{i-1})
 * where 
 *      r_{i+1/2} = (df/dx)_{i+1/2} / (df/dx)_{i-1/2}, 
 * and
 *      (df/dx)_{i+1/2} = (f_{i+1} - f_i)/(x_{i+1} - x_i).
 * Here, i+1/2 refers to the downstream cell face of the cell 
 * with centroid at x_i.
 * A few common schemes:
 *  First order upwind:  a = 0,   b = 0,
 *  Central difference:  a = 1,   b = 0,
 *  Second order upwind: a = 0,   b = 1,
 *  QUICK:               a = 3/4, b = 1/4.
 * 
 * Explicitly, the scheme is given by:
 * 
 * Positive advection coefficient:
 *   f_{i+1/2} =
 *       (1 + 0.5*dx_i * [ -a/dx_{i+1/2} + b/dx_{i-1/2} ]   ) * f_i
 *     + 0.5*a*dx_i/dx_{i+1/2} * f_{i+1}
 *     - 0.5*b*dx_i/dx_{i-1/2} * f_{i-1}
 * Negative advection coefficient:
 *   f_{i+1/2} =
 *       (1 + 0.5*dx_i * [ -a/dx_{i+1/2} + b/dx_{i-1/2} ]   ) * f_{i+1}
 *     + 0.5*a*dx_i/dx_{i+1/2} * f_i
 *     - 0.5*b*dx_i/dx_{i-1/2} * f_{i+2}
*/



/**
 * Interpolation coefficients with non-uniform QUICK (quadratic upwind interpolation):
 * a third-order upwind-biased (advection-sign-dependent) 3-point stencil.
 * delta includes the interpolation coefficients for the density q_{i-1/2}
 * evaluated on the cell interface between cells i-1 and i.
 */
void AdvectionDiffusionTerm::SetInterpolationCoefficientQUICK(real_t ***deltas, real_t **A, FVM::fluxGridType fluxGridType){
    len_t nr = this->AdvectionTerm::nr + (fluxGridType==FVM::FLUXGRIDTYPE_RADIAL);

    const real_t *dx, *dx_f;
    for(len_t ir=0; ir<nr; ir++){
        len_t n1,n2;
        // XXX: assume same momentum grid at all radii for radial flux grid
        if(fluxGridType==FVM::FLUXGRIDTYPE_RADIAL){
            n1 = this->AdvectionTerm::n1[0];
            n2 = this->AdvectionTerm::n2[0];
        } else {
            n2 = this->AdvectionTerm::n2[ir]+(fluxGridType==FVM::FLUXGRIDTYPE_P2);
            n1 = this->AdvectionTerm::n1[ir]+(fluxGridType==FVM::FLUXGRIDTYPE_P1);
        }
        MomentumGrid *mg = this->AdvectionTerm::grid->GetMomentumGrid(ir);
        int N;
        switch(fluxGridType){
            case FVM::FLUXGRIDTYPE_RADIAL:
                dx   = this->AdvectionTerm::grid->GetRadialGrid()->GetDr();
                dx_f = this->AdvectionTerm::grid->GetRadialGrid()->GetDr_f();
                N = nr;
                break;
            case FVM::FLUXGRIDTYPE_P1:
                dx   = mg->GetDp1();
                dx_f = mg->GetDp1_f();
                N = n1;
                break;
            case FVM::FLUXGRIDTYPE_P2:
                dx   = mg->GetDp2();
                dx_f = mg->GetDp2_f();
                N = n2;
                break;
            default:
                dx   = nullptr;
                dx_f = nullptr;
        }
        
        for(len_t i=0; i<n1; i++)
            for(len_t j=0; j<n2; j++){
                bool isFlowPositive = (A[ir][j*n1+i]>0);
                
                // sign of A: +1 for A>=0 and -1 for A<0 
                int sgn = isFlowPositive - !isFlowPositive;

                real_t dxU1, dxU2, dxD1, dx0_f, dxU1_f;
                int ind;
                switch(fluxGridType){
                    case FVM::FLUXGRIDTYPE_RADIAL:
                        ind = ir;
                        break;
                    case FVM::FLUXGRIDTYPE_P1:
                        ind = i;
                        break;
                    case FVM::FLUXGRIDTYPE_P2:
                        ind = j;
                        break;
                    default:
                        break;
                }

                int shiftU1 = -(1+sgn)/2;   // 0.5 step upstream
                int shiftU2 = -(1+3*sgn)/2; // 1.5 steps upstream
                int shiftD1 = (-1+sgn)/2;   // 0.5 step downstream


                /**
                 * All interpolation points are inside the grid when
                 * (3+sgn)/2 <= ind <= N - (5-sgn)/2.
                 * Otherwise, assume constant grid spacing
                 */
                if( (ind < 1 + (1+sgn)/2) || (ind > N-3 + (1+sgn)/2)) 
                    dxU1 = dxU2 = dxD1 = dx0_f = dxU1_f = 1;
                else {
                    dxU1 = dx[ind + shiftU1];
                    dxU2 = dx[ind + shiftU2];
                    dxD1 = dx[ind + shiftD1];
                    dx0_f  = dx_f[ind-1];
                    dxU1_f = dx_f[ind-1-sgn];
                }

                deltas[ir][j*n1+i][2+shiftU1] = dxD1*( dxU1 + 0.5*dxU2 ) / ( 2*dx0_f*dxU1_f );
                deltas[ir][j*n1+i][2+shiftU2] = - dxU1*dxD1 /( 4*dxU1_f*(dx0_f+dxU1_f) );
                deltas[ir][j*n1+i][2+shiftD1] = dxU1*( dxU1 + 0.5*dxU2 ) / ( 2*dx0_f*(dx0_f+dxU1_f) );


/*
                // Apply boundary conditions: 
                // deltas[k] is the interpolation weight for f_{i+k-2}
                // Mirrored boundary condition: f_{-|i|} = f_{|i|-1}, i>0
                // Dirichlet: f_{|i|} = 0 (default method in AdvectionTerm.set, ie no flux through boundary)

                if(ind+shiftU1 < 0){
                    if( bc_lower == AD_BC_MIRRORED )
                        deltas[ir][j*n1+i][1-2*ind - shiftU1] += deltas[ir][j*n1+i][2+shiftU1];
                    deltas[ir][j*n1+i][2+shiftU1] = 0;
                }
                if(ind+shiftU2 < 0){
                    if( bc_lower == AD_BC_MIRRORED )
                        deltas[ir][j*n1+i][1-2*ind - shiftU2] += deltas[ir][j*n1+i][2+shiftU2];
                    deltas[ir][j*n1+i][2+shiftU2] = 0;
                } 
                if(ind+shiftD1 < 0){
                    if( bc_lower == AD_BC_MIRRORED )
                        deltas[ir][j*n1+i][1-2*ind - shiftD1] += deltas[ir][j*n1+i][2+shiftD1];
                    deltas[ir][j*n1+i][2+shiftD1] = 0;
                }
                
                if(ind+shiftU1 > N-2) {
                    if( bc_upper == AD_BC_MIRRORED )
                        deltas[ir][j*n1+i][2*N-1-2*ind - shiftU1] += deltas[ir][j*n1+i][2+shiftU1];
                    deltas[ir][j*n1+i][2+shiftU1] = 0;
                }
                if(ind+shiftU2 > N-2) {
                    if( bc_upper == AD_BC_MIRRORED )
                        deltas[ir][j*n1+i][2*N-1-2*ind - shiftU2] += deltas[ir][j*n1+i][2+shiftU2];
                    deltas[ir][j*n1+i][2+shiftU2] = 0;
                }
                if(ind+shiftD1 > N-2) {
                    if( bc_upper == AD_BC_MIRRORED )
                        deltas[ir][j*n1+i][2*N-1-2*ind - shiftD1] += deltas[ir][j*n1+i][2+shiftD1];
                    deltas[ir][j*n1+i][2+shiftD1] = 0;
                }



mirror at 0: f_-|i| -> f_|i|-1, i>0
f_x = delta[x+2-ind] 
delta[k] ~ f_{ind+k-2} -> f_{1-k-ind} = delta[3-2*ind-k] 
delta[2+shift] ~ f_{ind+shift} -> f_{-ind-shift-1} ~ delta[1-2*ind-shift]

mirror at N: 
f_{N-2+x} -> f_{N-1-x}, x>0
f_i -> f_{2*N-3-i}
delta[2+shift] ~ f_{ind+shift} -> f_{2*N-3-ind-shift} ~ delta[2*N-1-2*ind-shift)]
*/
            }
    }
}

void AdvectionDiffusionTerm::ApplyInterpolationCoefficientBoundaryCondition(
    real_t ***deltas, adv_bc bc_lower,adv_bc bc_upper, FVM::fluxGridType fluxGridType){

    len_t nr = this->AdvectionTerm::nr + (fluxGridType==FVM::FLUXGRIDTYPE_RADIAL);    for(len_t ir=0; ir<nr; ir++){
        len_t n1,n2;
        // XXX: assume same momentum grid at all radii for radial flux grid
        if(fluxGridType==FVM::FLUXGRIDTYPE_RADIAL){
            n1 = this->AdvectionTerm::n1[0];
            n2 = this->AdvectionTerm::n2[0];
        } else {
            n2 = this->AdvectionTerm::n2[ir]+(fluxGridType==FVM::FLUXGRIDTYPE_P2);
            n1 = this->AdvectionTerm::n1[ir]+(fluxGridType==FVM::FLUXGRIDTYPE_P1);
        }

        for(len_t i=0; i<n1; i++)
            for(len_t j=0; j<n2; j++){
                len_t ind, N;
                switch(fluxGridType){
                    case FVM::FLUXGRIDTYPE_RADIAL:
                        ind = ir;
                        N = nr-1;
                        break;
                    case FVM::FLUXGRIDTYPE_P1:
                        ind = i;
                        N = n1-1;
                        break;
                    case FVM::FLUXGRIDTYPE_P2:
                        ind = j;
                        N = n2-1;
                        break;
                    default:
                        break;
                }


                if(ind==0){
                    if(bc_lower == AD_BC_MIRRORED){
                        deltas[ir][j*n1+i][3] += deltas[ir][j*n1+i][0];
                        deltas[ir][j*n1+i][2] += deltas[ir][j*n1+i][1];
                    }    
                    deltas[ir][j*n1+i][0] = 0;
                    deltas[ir][j*n1+i][1] = 0;
                } else if(ind==1) {
                    if(bc_lower == AD_BC_MIRRORED)
                        deltas[ir][j*n1+i][2] += deltas[ir][j*n1+i][0];
                    deltas[ir][j*n1+i][0] = 0;
                } 
                
                if (ind==N-1) {
                    if(bc_upper == AD_BC_MIRRORED){
                        deltas[ir][j*n1+i][1] += deltas[ir][j*n1+i][2];
                        deltas[ir][j*n1+i][0] += deltas[ir][j*n1+i][3];
                    }
                    deltas[ir][j*n1+i][2] = 0;
                    deltas[ir][j*n1+i][3] = 0;

                } else if (ind==N-2) {
                    if(bc_upper == AD_BC_MIRRORED)
                        deltas[ir][j*n1+i][2] += deltas[ir][j*n1+i][3];
                    deltas[ir][j*n1+i][3] = 0;
                }

            }
    }
}


/**
 * Set all interpolation coefficients to the same value.
 *
 * v: Value to assign to all interpolation coefficients.
 */
void AdvectionDiffusionTerm::SetInterpolationCoefficientValues(const real_t v) {
    for (len_t k = 0; k < this->AdvectionTerm::nr; k++) {
        const len_t n2 = this->AdvectionTerm::n2[k];
        const len_t n1 = this->AdvectionTerm::n1[k];

        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1; i++) {
                this->AdvectionTerm::deltar[k][j*n1 + i] = v;
                this->AdvectionTerm::delta1[k][j*n1 + i] = v;
                this->AdvectionTerm::delta2[k][j*n1 + i] = v;
            }
        }
    }

    const len_t nr = this->AdvectionTerm::nr;

    for (len_t ir = 0; ir < nr; ir++) {
        const len_t n2 = this->AdvectionTerm::n2[ir];
        const len_t n1 = this->AdvectionTerm::n1[ir];
        // reset coefficients
        for(len_t i=0; i<(n1+1)*n2; i++)
            for(len_t k=0; k<2*stencil_order; k++)
                delta1s[ir][i][k] = 0;
        for(len_t i=0; i<n1*(n2+1); i++)
            for(len_t k=0; k<2*stencil_order; k++)
                delta2s[ir][i][k] = 0;

        for (len_t j = 0; j < n2; j++){
            for (len_t i = 1; i < n1; i++) {
                delta1s[ir][j*(n1+1)+i][stencil_order-1] = 1-v;
                delta1s[ir][j*(n1+1)+i][stencil_order] = v;
            }

            // Set dirichlet on upper, neumann on lower (internal) boundary
            delta1s[ir][j*(n1+1)][stencil_order] = 1;
        }
        for (len_t i = 0; i < n1; i++) {
            for (len_t j = 1; j < n2; j++){
                delta2s[ir][j*n1+i][stencil_order-1] = 1-v;
                delta2s[ir][j*n1+i][stencil_order] = v;
            }        
            // Set neumann on both boundaries (assumed to be internal)
            delta2s[ir][i][stencil_order] = 1;
            delta2s[ir][n2*n1+i][stencil_order] = 1;
        }
    }


    const len_t n2 = this->AdvectionTerm::n2[0]; // XXX: assume same grid at all radii
    const len_t n1 = this->AdvectionTerm::n1[0];
    // reset coefficients
    for(len_t ir=0; ir<nr+1; ir++)
        for(len_t i=0; i<n1*n2; i++)
            for(len_t k=0;k<2*stencil_order;k++)
                deltars[ir][i][k] = 0;

    for (len_t ir = 1; ir < nr; ir++) {
        for (len_t i = 0; i < n1*n2; i++) {
            deltars[ir][i][stencil_order-1] = 1-v;
            deltars[ir][i][stencil_order] = v;
        }
    }

    // set dirichlet on upper, neumann on lower (internal)
    for (len_t i = 0; i < n1*n2; i++) 
        deltars[0][i][stencil_order] = 1;
}






/**
 * Returns the delta coefficient given a Peclet number Pe,
 * going smoothly between a central difference scheme in diffusion
 * dominated cases to first-order upwind difference in advection
 * dominated cases (the later characterised by |Pe|>>1) 
 * Solutions are quite sensitive to Pe_threshold, and
 * there is no recipe for how to choose it in general. 
 * Somewhere between 2 and 100 seems good in some cases.
 */
real_t UpwindDelta(real_t A, real_t D, real_t dx){
    real_t Pe_threshold = 10;
    real_t delta_central = 0.5;        
    real_t delta_upwind = (A<0)  - (A>0);

    if(D==0) // upwind for pure advection
        return delta_upwind;

    real_t Pe = dx*A/D;

    real_t limiter;
    real_t x = abs(Pe)/Pe_threshold;
    if(x<=1)
        limiter = 1;
    else
        limiter = 1 / (1 + (x-1)*(x-1)/(x*sqrt(x)) );
    return limiter*delta_central + (1-limiter)*delta_upwind;
}

/** 
 * Set interpolation coefficients dynamically based on the Peclet number
 * of the advection and diffusion coefficients currently stored in memory.
 * Upwind difference for |Pe|>Pe_threshold and central otherwise.
 * Pe_threshold should be at least two, below which central difference 
 * is known to always be stable.
 */
void AdvectionDiffusionTerm::SetInterpolationCoefficientValuesHybridUpwind() {
    real_t fr, f1, f2, drr, d11, d22;
    real_t DeltaR, DeltaP1, DeltaP2;
    for (len_t ir = 0; ir < this->AdvectionTerm::nr; ir++) {
        const len_t n2 = this->AdvectionTerm::n2[ir];
        const len_t n1 = this->AdvectionTerm::n1[ir];

        for (len_t j = 0; j < n2; j++) {
            for (len_t i = 0; i < n1; i++) {
                fr  = this->AdvectionTerm::fr[ir][j*n1+i] ;
                f1  = this->AdvectionTerm::f1[ir][j*(n1+1)+i] ;
                f2  = this->AdvectionTerm::f2[ir][j*n1+i] ;
                drr = this->DiffusionTerm::drr[ir][j*n1+i];
                d11 = this->DiffusionTerm::d11[ir][j*(n1+1)+i];
                d22 = this->DiffusionTerm::d22[ir][j*n1+i];
                DeltaR  = this->AdvectionTerm::grid->GetRadialGrid()->GetDr(ir);
                DeltaP1 = this->AdvectionTerm::grid->GetMomentumGrid(ir)->GetDp1(i);
                DeltaP2 = this->AdvectionTerm::grid->GetMomentumGrid(ir)->GetDp2(j);

                this->AdvectionTerm::deltar[ir][j*n1 + i] = UpwindDelta(fr,drr,DeltaR);
                this->AdvectionTerm::delta1[ir][j*n1 + i] = UpwindDelta(f1,d11,DeltaP1);
                this->AdvectionTerm::delta2[ir][j*n1 + i] = UpwindDelta(f2,d22,DeltaP2);
            }
        }
    }


}

/**
 * Save the advection and diffusion coefficients of this
 * object to the specified file.
 */
void AdvectionDiffusionTerm::SaveCoefficientsSFile(const string& filename) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_WRITE);
    this->SaveCoefficientsSFile(sf);
    sf->Close();
}
void AdvectionDiffusionTerm::SaveCoefficientsSFile(SFile *sf) {
    this->AdvectionTerm::SaveCoefficientsSFile(sf);
    this->DiffusionTerm::SaveCoefficientsSFile(sf);
}

