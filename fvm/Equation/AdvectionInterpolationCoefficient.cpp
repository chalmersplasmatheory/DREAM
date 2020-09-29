/**
 * Implementation of a class which evaluates, stores and provides
 * interpolation coefficients for advection terms.
 */
#include "FVM/Equation/AdvectionInterpolationCoefficient.hpp"
#include <algorithm>
#include <limits>

using namespace DREAM::FVM;
using namespace std;



/**
 * Constructor.
 */
AdvectionInterpolationCoefficient::AdvectionInterpolationCoefficient(Grid*g, fluxGridType fluxGridType, 
        adv_bc bc_l, adv_bc bc_u){
    this->grid = g;
    this->fgType = fluxGridType;
    this->bc_lower = bc_l;
    this->bc_upper = bc_u;

    this->delta_prev = new real_t[2*stencil_width];
    delta_prev[0] = -1; // indicator that it is uninitialised
    for(len_t k=1; k<2*stencil_width;k++)
        delta_prev[k] = 0;
}

/**
 * Destructor.
 */
AdvectionInterpolationCoefficient::~AdvectionInterpolationCoefficient(){
    Deallocate();
    delete [] delta_prev;
}


/**
 * To be called when grids are rebuilt; allocates and initializes interpolation coefficients
 */
bool AdvectionInterpolationCoefficient::GridRebuilt(){
    Deallocate();
    this->nr = grid->GetNr() + (fgType==FLUXGRIDTYPE_RADIAL);
    n1 = new len_t[nr];
    n2 = new len_t[nr];
    deltas = new real_t**[nr];
    deltas_jac = new real_t**[nr];
    for(len_t ir=0; ir<nr; ir++){
        // XXX: If radial flux grid, assume same momentum grid at all radii
        n1[ir] = grid->GetNp1(ir*(fgType!=FLUXGRIDTYPE_RADIAL)) + (fgType==FLUXGRIDTYPE_P1);
        n2[ir] = grid->GetNp2(ir*(fgType!=FLUXGRIDTYPE_RADIAL)) + (fgType==FLUXGRIDTYPE_P2);
        deltas[ir] = new real_t*[n1[ir]*n2[ir]];
        deltas_jac[ir] = new real_t*[n1[ir]*n2[ir]];
        for(len_t i=0; i<n1[ir]*n2[ir]; i++){
            deltas[ir][i] = new real_t[2*stencil_width];
            deltas_jac[ir][i] = new real_t[2*stencil_width];
        }
        
    }
    ResetCoefficient();
    return true;
}

/**
 * Sets interpolation coefficients to 0.
 */
void AdvectionInterpolationCoefficient::ResetCoefficient(){
    for(len_t ir=0; ir<nr; ir++)
        for(len_t i=0; i<n1[ir]*n2[ir]; i++)
            for(len_t k=0; k<2*stencil_width; k++){
                deltas[ir][i][k] = 0;
                deltas_jac[ir][i][k] = 0;
            }
    hasNonTrivialJacobian = false;
}


/**
 * Set the interpolation coefficients delta based on interpolation method adv_i
 */
void AdvectionInterpolationCoefficient::SetCoefficient(real_t **A, real_t **/*D*/, UnknownQuantityHandler *unknowns, adv_interpolation adv_i, real_t damping_factor){
    if(!hasBeenInitialized)
        hasBeenInitialized = GridRebuilt();
    else
        ResetCoefficient();
    SetNNZ(adv_i);
    const real_t *x = nullptr, *x_f = nullptr;
    int_t N;
    for(len_t ir=0; ir<nr; ir++){
        switch(fgType){
            case FLUXGRIDTYPE_RADIAL:
                x   = grid->GetRadialGrid()->GetR();
                x_f = grid->GetRadialGrid()->GetR_f();
                break;
            case FLUXGRIDTYPE_P1:
                x   = grid->GetMomentumGrid(ir)->GetP1();
                x_f = grid->GetMomentumGrid(ir)->GetP1_f();
                break;
            case FLUXGRIDTYPE_P2:
                x   = grid->GetMomentumGrid(ir)->GetP2();
                x_f = grid->GetMomentumGrid(ir)->GetP2_f();
            default:
                break;
        }
        
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                len_t pind = j*n1[ir]+i;

                bool isFlowPositive = (A[ir][pind]>0);
                
                // sign of A: +1 for A>0 and -1 for A<=0 
                int_t sgn = isFlowPositive - !isFlowPositive;
                shiftU1 = -(1+sgn)/2;   // 0.5 step upstream
                shiftU2 = -(1+3*sgn)/2; // 1.5 steps upstream
                shiftD1 = (-1+sgn)/2;   // 0.5 step downstream
                
                int_t ind = GetIndex(ir,i,j,&N);
                xf  = x_f[ind];
                x_0 = x_f[0];
                xN  = x_f[N];

                adv_interpolation method = adv_i;
                // When 1 or 2 grid points are used, use central difference scheme
                if(N<3)
                    method = AD_INTERP_CENTRED;
                // on the initial rebuild, if using flux limiters, use a robust method 
                // to obtain a good starting point to iterate the flux limiter from
                else if(isFirstRebuild && IsFluxLimiterMethod(adv_i))
                    method = AD_INTERP_UPWIND;


                real_t alpha;
                switch(method){
                    case AD_INTERP_CENTRED: {
                        alpha = 0.0;
                        SetFirstOrderCoefficient(ind,N,x,alpha,deltas[ir][pind]);
                        break;
                    } case AD_INTERP_UPWIND: {
                        alpha = 0.5;
                        SetFirstOrderCoefficient(ind,N,x,alpha,deltas[ir][pind]);
                        break;
                    } case AD_INTERP_UPWIND_2ND_ORDER: {
                        // 2nd order upwind
                        alpha = -(1.0/8.0)*(GetXi(x,ind+shiftU2,N) - xf)/(GetXi(x,ind+shiftD1,N) - xf);
                        SetSecondOrderCoefficient(ind,N,x,alpha,deltas[ir][pind]);
                        break;
                    } case AD_INTERP_DOWNWIND: {
                        alpha = 0.5*(GetXi(x,ind+shiftD1,N) - xf)/(GetXi(x,ind+shiftU1,N) - xf);
                        SetFirstOrderCoefficient(ind,N,x,alpha,deltas[ir][pind]);
                        break;
                    } case AD_INTERP_QUICK: {
                        alpha = 0.0;
                        SetSecondOrderCoefficient(ind,N,x,alpha,deltas[ir][pind]);
                        break;
                    } case AD_INTERP_SMART: {
                        // Sets interpolation coefficients using the flux limited
                        // SMART method

                        std::function<real_t(int_t)> yFunc = GetYFunc(ir,i,j,unknowns);
                        real_t r = GetFluxLimiterR(ind,N,yFunc,x);
                        
                        real_t kappa = 0.5;
                        real_t M = 4;
                        alpha = 0;
                        SetGPLKScheme(ind, N, x, r, alpha, kappa, M, damping_factor, deltas[ir][pind]);

                        break;
                    } case AD_INTERP_MUSCL: {
                        // Sets interpolation coefficients using the flux limited
                        // MUSCL method

                        std::function<real_t(int_t)> yFunc = GetYFunc(ir,i,j,unknowns);
                        real_t r = GetFluxLimiterR(ind,N,yFunc,x);
                        
                        real_t kappa = 0;
                        real_t M = 2;
                        alpha = 0;
                        SetGPLKScheme(ind, N, x, r, alpha, kappa, M, damping_factor, deltas[ir][pind]);

                        break;
                    } case AD_INTERP_OSPRE: {
                        // Sets interpolation coefficients using the continuous
                        // flux limited OSPRE method

                        std::function<real_t(int_t)> yFunc = GetYFunc(ir,i,j,unknowns);
                        real_t r = GetFluxLimiterR(ind,N,yFunc,x);
                        real_t psi = 1.5*r*(r+1)/(r*r+r+1);
//                        real_t psiPrime = 1.5*(1+2*r)/((r*r+r+1)*(r*r+r+1));
//                        real_t a = psi - r*psiPrime;
//                        real_t b = psiPrime;
//                        SetLinearFluxLimitedCoefficient(ind,N,x,a,b,deltas[ir][pind]);
                        SetFluxLimitedCoefficient(ind,N,x,psi,deltas[ir][pind]);
//                        SetFluxLimitedCoefficient(ind,N,x,psi,deltas_jac[ir][pind],r,psiPrime);
//                        hasNonTrivialJacobian = true;

                        break;
                    } case AD_INTERP_TCDF: {
                        // Sets interpolation coefficients using the continuous
                        // flux limited TCDF method, described in
                        //      D Zhang et al., J Comp Phys 302, 114 (2015).
                        // The flux limiter is defined piecewise but is C1 in order
                        // to ensure good convergence properties. Has a large overlap 
                        // with QUICK (in the interval 0.5 < r < 2.0)

                        std::function<real_t(int_t)> yFunc = GetYFunc(ir,i,j,unknowns);
                        real_t r = GetFluxLimiterR(ind,N,yFunc,x);
                        real_t psi, psiPrime;
                        if(r<0){
                            psi = r*(1+r)/(1+r*r);
                            psiPrime = (1+2*r-r*r)/((1+r*r)*(1+r*r));
                        } else if(r<0.5){
                            psi = r*r*r - 2*r*r +2*r;
                            psiPrime = 3*r*r - 4*r + 2;
                        } else if(r<2){
                            psi = 0.25 + 0.75*r;
                            psiPrime = 0.75;
                        } else {
                            psi = (2*r*r - 2*r - 2.25) / (r*r - r -1);
                            psiPrime = 0.25*(2*r-1) / ((r*r - r -1)*(r*r - r -1)); 
                        }
                        SetFluxLimitedCoefficient(ind,N,x,psi,deltas[ir][pind]);
                        SetFluxLimitedCoefficient(ind,N,x,psi,deltas_jac[ir][pind],r,psiPrime);
                        hasNonTrivialJacobian = true;
/*
                        real_t a = psi - r*psiPrime;
                        real_t b = psiPrime;
                        SetLinearFluxLimitedCoefficient(ind,N,x,a,b,deltas[ir][pind]);
*/
                        break;
                    } default: {
                        throw FVMException("Invalid interpolation method: not yet supported.");
                        break;
                    }
                }

                if(!hasNonTrivialJacobian)
                    for(len_t k=0;k<2*stencil_width; k++)
                        deltas_jac[ir][pind][k] = deltas[ir][pind][k];

                // set nearly zero interpolation coefficients to identically zero
                real_t eps = std::numeric_limits<real_t>::epsilon();
                real_t threshold_eps = 1e6;
                for(len_t k=0; k<2*stencil_width; k++){
                    if(abs(deltas[ir][pind][k]) < eps*threshold_eps)
                        deltas[ir][pind][k] = 0.0;
                    if(abs(deltas_jac[ir][pind][k]) < eps*threshold_eps)
                        deltas_jac[ir][pind][k] = 0.0;
                }
                isFirstRebuild = false;
            }
    }
    ApplyBoundaryCondition();
}

/**
 * Sets flux limiter according to the "generalized piecewise linear kappa-scheme" with
 * parameters kappa, alpha, M, according to the nomenclature in 
 *      N P Waterson, H Deconinck, JCP 224 (2007).
 * Requirements for boundedness: -1<=kappa<=1, -1<=alpha<=0, 1<=M.
 *  MUSCL: kappa=0;       M=2; alpha=0
 *  SMART: kappa=0.5;     M=4; alpha=0
 *  KOREN: kappa=1.0/3.0; M=2; alpha=0. (not implemented)
 * Support is added for extending the range of the kappa scheme due to finite Peclet numbers
 * (since diffusion helps stabilize), following 
 *      H Smaoui et al, Int. J. Comp. Meth. Eng. Sci. Mech. 9, 180 (2008)
 * but testing shows that overall accuracy does not improve over the regular (PeInv=0) schemes
 */
void AdvectionInterpolationCoefficient::SetGPLKScheme(int_t ind, int_t N, const real_t *x, real_t r, real_t alpha, real_t kappa, real_t M, real_t damping, real_t *&deltas){
    real_t a0 = 0.5*(1-kappa);
    real_t b0 = 0.5*(1+kappa);
    real_t b1 = 2.0 + alpha;
    real_t a2 = M;
    // real_t B3 = 2*PeInv;
    if(delta_prev[0] == -1){ // initialize with the target kappa scheme
        delta_prev[0] = 0;
//        SetFluxLimitedCoefficient(ind,N,x,a0+b0*r,delta_prev);
        SetLinearFluxLimitedCoefficient(ind,N,x,a0,b0,delta_prev);
    } 
    real_t a,b;
    if(r<=0) {
        a = 0;
        b = 0;
    } else if ( (b1*r < a2) && (b1*r < a0+b0*r) ) {
        a = 0;
        b = b1;
    } else if (a0+b0*r <= a2) {
        a = a0;
        b = b0;
    } else {
        a = a2;
        b = 0;
    }
//    SetFluxLimitedCoefficient(ind,N,x,a+b*r,deltas);
    SetLinearFluxLimitedCoefficient(ind,N,x,a,b,deltas);
    for(len_t k=0; k<2*stencil_width; k++){
        deltas[k] = delta_prev[k] + damping * (deltas[k] - delta_prev[k]); 
        delta_prev[k] = deltas[k];
    }
    
}


/**
 * Sets second order interpolation schemes:
 *      alpha = 0: QUICK (3rd order)
 *      alpha = -1/8: Central difference
 *      alpha = -b/(8*c): Second order upwind
 */                        
void AdvectionInterpolationCoefficient::SetSecondOrderCoefficient(int_t ind, int_t N, const real_t *x, real_t alpha, real_t *&deltas){
    real_t a = GetXi(x,ind+shiftU1,N) - xf;
    real_t b = GetXi(x,ind+shiftU2,N) - xf;
    real_t c = GetXi(x,ind+shiftD1,N) - xf;

    deltas[2+shiftU1] = c*(b+8*alpha*a)/( (a-b)*(a-c) );
    deltas[2+shiftU2] = -(1+8*alpha)*a*c/( (a-b)*(b-c) );
    deltas[2+shiftD1] = a*(b+8*alpha*c)/( (a-c)*(b-c) );
}

/**
 * Sets first order interpolation schemes:
 *      alpha = 0: Central difference (2nd order)
 *      alpha = 0.5: First-order upwind
 *      alpha = 0.5*b/a: First-order downwind
 */
void AdvectionInterpolationCoefficient::SetFirstOrderCoefficient(int_t ind, int_t N, const real_t *x, real_t alpha, real_t *&deltas, real_t scaleFactor){
    real_t a = GetXi(x,ind+shiftU1,N) - xf;
    real_t b = GetXi(x,ind+shiftD1,N) - xf;

    if( b==a ) { 
        // in this case the grid is probably empty and the coefficient 
        // will not be used; defaulting it to first-order upwind
        deltas[2+shiftU1] = 1;
    } else {
        deltas[2+shiftU1] = scaleFactor* (b-2*alpha*a)/(b-a);
        deltas[2+shiftD1] = -scaleFactor* (1-2*alpha)*a/(b-a);
    }
}

/**
 * Sets non-linear flux limited interpolation schemes:
 *      y_{i-1/2} = y_{i-1} + (x_{i-1/2} - x_{i-1})*psi(r_{i-1/2})*(y_{i-1} - y_{i-2})/(x_{i-1} - x_{i-2}),
 *                = y_{i-1} + 0.5*psi(r_{i-1/2})*dx_{i-3/2} * y'_{i-3/2} 
 * where
 *      r_{i-1/2} = y'_{i-1/2} / y'_{i-3/2},
 * and
 *      y'_{i-1/2} = (y_i - y_{i-1})/(x_i - x_{i-1})
 *      y'_{i-3/2} = (y_{i-1} - y_{i-2})/(x_{i-1} - x_{i-2})
 */                        
void AdvectionInterpolationCoefficient::SetFluxLimitedCoefficient(int_t ind, int_t N, const real_t *x, real_t psi, real_t *&deltas, real_t r, real_t psiPrime){
    real_t x1 = GetXi(x,ind+shiftU1,N);
    real_t x2 = GetXi(x,ind+shiftU2,N);
    real_t k = (xf-x1)/(x1-x2); // = -/+ 0.5 for uniform
    deltas[2+shiftU1] = 1 + k*psi;
    deltas[2+shiftU2] = -k*psi;

    if(psiPrime==0)
        return;

    // if psiPrime is provided, add jacobian
    real_t x0 = GetXi(x,ind+shiftD1,N);    
    real_t l = (x1-x2)/(x0-x1); // = +/- 1 for uniform

    deltas[2+shiftD1] = k*psiPrime*l;
    deltas[2+shiftU1] += -k*psiPrime*(r+l);
    deltas[2+shiftU2] += k*psiPrime*r;
}


/**
 * Sets non-linear flux limited interpolation scheme, using a linear flux limiter function of the form
 *      psi(r) = a_psi + b_psi * r.
 * Unlike the general flux limiter coefficient, this method will have the "full" contribution to the jacobian
 */                        
void AdvectionInterpolationCoefficient::SetLinearFluxLimitedCoefficient(int_t ind, int_t N, const real_t *x, real_t a_psi, real_t b_psi, real_t *&deltas){
    real_t x0 = GetXi(x,ind+shiftD1,N);
    real_t x1 = GetXi(x,ind+shiftU1,N);
    real_t x2 = GetXi(x,ind+shiftU2,N);

    real_t dx0 = xf - x1;

    deltas[2+shiftU1] = 1 + a_psi * dx0/(x1-x2);
    deltas[2+shiftU2] = -a_psi * dx0/(x1-x2);
    deltas[2+shiftD1] = b_psi * dx0/(x0-x1);
    deltas[2+shiftU1] -= b_psi * dx0/(x0-x1);
}


/**
 * Apply default boundary conditions. 
 *      Mirrored: assumes that the solution is symmetric around the grid boundary
 *      Dirichlet: Zero flux through the boundary (equivalent to y=0 on boundary)
 */
void AdvectionInterpolationCoefficient::ApplyBoundaryCondition(){
    for(len_t ir=0; ir<nr; ir++)
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                int_t N;
                int_t ind = GetIndex(ir,i,j,&N);
                len_t pind = j*n1[ir]+i;
                len_t k_max = 2*stencil_width-1;
                if(bc_lower == AD_BC_MIRRORED){
                    for(len_t k=0; k+ind<stencil_width; k++){
                        deltas[ir][pind][k_max-2*ind-k] += deltas[ir][pind][k];
                        deltas[ir][pind][k] = 0;
                        deltas_jac[ir][pind][k_max-2*ind-k] += deltas_jac[ir][pind][k];
                        deltas_jac[ir][pind][k] = 0;
                    }
                } else if(bc_lower == AD_BC_DIRICHLET)
                    if(ind==0)
                        for(len_t k=0; k<2*stencil_width; k++){
                            deltas[ir][pind][k] = 0;
                            deltas_jac[ir][pind][k] = 0;
                        }
 
                if(bc_upper == AD_BC_MIRRORED){
                    for(len_t k=N+stencil_width-ind; k<=k_max; k++){
                        deltas[ir][pind][k_max+2*(N-ind)-k] += deltas[ir][pind][k];
                        deltas[ir][pind][k] = 0;
                        deltas_jac[ir][pind][k_max+2*(N-ind)-k] += deltas_jac[ir][pind][k];
                        deltas_jac[ir][pind][k] = 0;
                    }
                } else if(bc_upper == AD_BC_DIRICHLET)
                    if(ind==N)
                        for(len_t k=0; k<2*stencil_width; k++){
                            deltas[ir][pind][k] = 0;
                            deltas_jac[ir][pind][k] = 0;
                        }
            }
}


/**
 * Lower limit on summation index for interpolation in AdvectionTerm.set
 * (goes from i-stencil_order unless this takes you outside the grid)
 */
len_t AdvectionInterpolationCoefficient::GetKmin(len_t i, len_t *n){
    if( stencil_width > i ){
        *n = stencil_width-i;
        return 0;
    }
    else {
        *n = 0;
        return i-stencil_width;
    }
}
/**
 * Upper limit on summation index for interpolation in AdvectionTerm.set
 * (goes to i+stencil_order-1 unless this takes you outside the grid)
 */
len_t AdvectionInterpolationCoefficient::GetKmax(len_t i, len_t N){
    if( N > (i+stencil_width) )
        return i+stencil_width-1;
    else
        return N-1;
}

/**
 * Returns a lambda function y(ind) that returns the unknown quantity y 
 * evaluated at index ind (with the other indices given by ir, i and/or j)
 */
std::function<real_t(int_t)> AdvectionInterpolationCoefficient::GetYFunc(len_t ir, len_t i, len_t j, FVM::UnknownQuantityHandler *unknowns){
    len_t offset=0;
    if(fgType==FVM::FLUXGRIDTYPE_RADIAL){
        return [this,unknowns,i,j](int_t ind)
        {
            len_t offset = 0;
            for(int_t k=0; k<ind;k++)
                offset += n1[k]*n2[k];
            return unknowns->GetUnknownData(id_unknown)[offset+j*n1[ind]+i]; 
        };
    } else if (fgType==FVM::FLUXGRIDTYPE_P1){
        for(len_t k=0; k<ir;k++)
            offset+=(n1[k]-1)*n2[k];
        return [this,unknowns,ir,j,offset](int_t ind)
            {return unknowns->GetUnknownData(id_unknown)[offset+j*(n1[ir]-1)+ind];};
    } else {
        for(len_t k=0; k<ir;k++)
                offset+=n1[k]*(n2[k]-1);        
        return [this,unknowns,ir,i,offset](int_t ind)
            {return unknowns->GetUnknownData(id_unknown)[offset+ind*n1[ir]+i];};
    }
}

/**
 * Sets the nnz parameter based on interpolation scheme.
 */
void AdvectionInterpolationCoefficient::SetNNZ(adv_interpolation adv_i){
    bool isFirstOrder = ( (adv_i==AD_INTERP_CENTRED) || (adv_i==AD_INTERP_DOWNWIND) 
                       || (adv_i==AD_INTERP_UPWIND) );
    if(isFirstOrder)
        nnzPerRow = 6*1+1; // = 7
    else
        nnzPerRow = 6*stencil_width+1; // = 13
}

/**
 * Deallocator
 */
void AdvectionInterpolationCoefficient::Deallocate(){
    if(deltas != nullptr){
        for(len_t ir=0; ir<nr; ir++){
            for(len_t i=0; i<n1[ir]*n2[ir]; i++)
                delete [] deltas[ir][i];
            delete [] deltas[ir];
        }
        delete [] deltas;
    }
    if(deltas_jac != nullptr){
        for(len_t ir=0; ir<nr; ir++){
            for(len_t i=0; i<n1[ir]*n2[ir]; i++)
                delete [] deltas_jac[ir][i];
            delete [] deltas_jac[ir];
        }
        delete [] deltas_jac;
    }

    if(n1!=nullptr){
        delete [] n1;
        delete [] n2;
    }

}


/**
 * Returns x[i] unless i falls outside the grid,
 * in which case we extrapolate assuming the grid 
 * is mirrored around the boundaries, in the sense
 * ie x[-|i|] - x_f[0] = x_f[0] - x[|i|-1] and
 * x[N-1 + |i|] - x_f[N] = x_f[N] - x[N-|i|] 
 * for |i|>0. 
 * x contains the grid point and has length N,
 * and xN is the upper boundary xN=x_f[N].
 */
real_t AdvectionInterpolationCoefficient::GetXi(const real_t *x, int_t i, int_t N){
    if(i<0)
        return 2*x_0 - x[-i-1];
    else if(i>N-1)
        return 2*xN - x[2*N-i-1];
    else
        return x[i];
}

/**
 * Returns y[i] unless i falls outside the grid,
 * in which case we either return y at the mirrored 
 * grid point or 0, depending on boundary condition.
 */
real_t AdvectionInterpolationCoefficient::GetYi(int_t i, int_t N, std::function<real_t(int_t)> y){
    bool isLoMirrored  = (bc_lower==AdvectionInterpolationCoefficient::AD_BC_MIRRORED);
    bool isLoDirichlet = (bc_lower==AdvectionInterpolationCoefficient::AD_BC_DIRICHLET);
    bool isUpMirrored  = (bc_upper==AdvectionInterpolationCoefficient::AD_BC_MIRRORED);
    bool isUpDirichlet = (bc_upper==AdvectionInterpolationCoefficient::AD_BC_DIRICHLET);

    if(i<0){
        if(isLoMirrored)
            return y(-i-1);
        else if(isLoDirichlet)
            return 0; 
        else 
            throw FVMException("The provided advection interpolation coefficent lower boundary condition is not supported by the SMART scheme.");
    } else if(i>N-1){
        if(isUpMirrored)
            return y(2*N-1-i);
        else if(isUpDirichlet)
            return 0; 
        else 
            throw FVMException("The provided advection interpolation coefficent upper boundary condition is not supported by the SMART scheme.");
    } else
        return y(i);
}

/**
 * Returns the flux-limiter parameter
 *      r_{i-1/2} = y'_{i-1/2} / y'_{i-3/2},
 * where
 *      y'_{i-1/2} = (y_i - y_{i-1})/(x_i - x_{i-1})
 *      y'_{i-3/2} = (y_{i-1} - y_{i-2})/(x_{i-1} - x_{i-2})

 */
real_t AdvectionInterpolationCoefficient::GetFluxLimiterR(int_t ind, int_t N, std::function<real_t(int_t)> y, const real_t *x){
    int_t i0 = ind+shiftD1;
    int_t i1 = ind+shiftU1;
    int_t i2 = ind+shiftU2;

    real_t x0 = GetXi(x, i0, N);
    real_t x1 = GetXi(x, i1, N);
    real_t x2 = GetXi(x, i2, N);
    
    real_t y0 = GetYi(i0, N, y);
    real_t y1 = GetYi(i1, N, y);
    real_t y2 = GetYi(i2, N, y);

    real_t dy0 = (y0-y1)/(x0-x1);
    real_t dy1 = (y1-y2)/(x1-x2);

    if(dy1==0) // return "essentially inifinity" with the sign of dy0
        return 1e5*( (dy0>0) - (dy0<0));

    real_t r = dy0/dy1;

    if(abs(r)>1e5)
        r = 1e5 * ( (r>0) - (r<0) ); 
    return r;
}
