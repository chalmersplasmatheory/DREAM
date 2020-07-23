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
    for(len_t ir=0; ir<nr; ir++){
        // XXX: If radial flux grid, assume same momentum grid at all radii
        n1[ir] = grid->GetNp1(ir*(fgType!=FLUXGRIDTYPE_RADIAL)) + (fgType==FLUXGRIDTYPE_P1);
        n2[ir] = grid->GetNp2(ir*(fgType!=FLUXGRIDTYPE_RADIAL)) + (fgType==FLUXGRIDTYPE_P2);
        deltas[ir] = new real_t*[n1[ir]*n2[ir]];
        for(len_t i=0; i<n1[ir]*n2[ir]; i++)
            deltas[ir][i] = new real_t[2*stencil_width];
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
            for(len_t k=0; k<2*stencil_width; k++)
                deltas[ir][i][k] = 0;
}





/**
 * Set the interpolation coefficients delta based on interpolation method adv_i
 */
void AdvectionInterpolationCoefficient::SetCoefficient(real_t **A, UnknownQuantityHandler *unknowns, adv_interpolation adv_i, real_t damping_factor){
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
                
                real_t alpha=0.0;
                // When 1 or 2 grid points are used, use central difference scheme
                if(N<3){
                    SetFirstOrderCoefficient(ind,N,x,alpha,deltas[ir][pind]);
                    continue;
                }
                switch(adv_i){
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
                        if(delta_prev[0] == -1){ // initialize with QUICK
                            delta_prev[0] = 0;
                            SetSecondOrderCoefficient(ind,N,x,0.0,delta_prev);
                        }                        
                        std::function<real_t(int_t)> yFunc = GetYFunc(ir,i,j,unknowns);
                        real_t phi = GetPhiHatNV(ind,N,yFunc);
                        if( (phi>1) || (phi<0) ){
                            // y_{i-1/2} = y_{i-1}: Upwind
                            alpha = 0.5;
                            SetFirstOrderCoefficient(ind,N,x,alpha,deltas[ir][pind]);
                        } else if (phi < 1.0/6.0) {
                            // y_{i-1/2} = 3*y_{i-1}: Amplified upwind
                            alpha = 0.5;
                            real_t scaleFactor = 3.0;
                            SetFirstOrderCoefficient(ind,N,x,alpha,deltas[ir][pind],scaleFactor);
                        } else if (phi > 5.0/6.0) {
                            // y_{i-1/2} = y_{i}: Downwind
                            alpha = 0.5*(GetXi(x,ind+shiftD1,N) - xf)/(GetXi(x,ind+shiftU1,N) - xf);
                            SetFirstOrderCoefficient(ind,N,x,alpha,deltas[ir][pind]);
                        } else {
                            // QUICK
                            alpha = 0.0;
                            SetSecondOrderCoefficient(ind,N,x,alpha,deltas[ir][pind]);
                        }

                        // APPLY DAMPING 
                        for(len_t k=0; k<2*stencil_width; k++){
                            deltas[ir][pind][k] = delta_prev[k] + damping_factor * (deltas[ir][pind][k] - delta_prev[k]); 
                            delta_prev[k] = deltas[ir][pind][k];
                        }

                        break;
                    } case AD_INTERP_MUSCL: {
                        // Sets interpolation coefficients using the flux limited
                        // MUSCL method

                        std::function<real_t(int_t)> yFunc = GetYFunc(ir,i,j,unknowns);
                        real_t r = GetFluxLimiterR(ind,N,yFunc,x);
                        if(delta_prev[0] == -1){ // initialize with Fromm
                            delta_prev[0] = 0;
                            SetLinearFluxLimitedCoefficient(ind,N,x,0.5,0.5,delta_prev);
                        }
                        
                        real_t a,b;
                        if(r<=0){
                            a = 0.0;
                            b = 0.0;
                        } else if (r>=3){
                            a = 2.0;
                            b = 0.0;
                        } else if (r>=1.0/3.0){
                            a = 0.5; 
                            b = 0.5;
                        } else {
                            a = 0.0;
                            b = 2.0;
                        }
                        SetLinearFluxLimitedCoefficient(ind,N,x,a,b,deltas[ir][pind]);
                        for(len_t k=0; k<2*stencil_width; k++){
                            deltas[ir][pind][k] = delta_prev[k] + damping_factor * (deltas[ir][pind][k] - delta_prev[k]); 
                            delta_prev[k] = deltas[ir][pind][k];
                        }

                        break;
                    } default: {
                        throw FVMException("Invalid interpolation method: not yet supported.");
                        break;
                    }

                    // set nearly zero interpolation coefficients to identically zero
                    real_t eps = std::numeric_limits<real_t>::epsilon();
                    for(len_t k=0; k<2*stencil_width; k++)
                        if(abs(deltas[ir][pind][k]) < 1e2*eps)
                            deltas[ir][pind][k] = 0.0;


                }
            }
    }
    ApplyBoundaryCondition();
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

    deltas[2+shiftU1] = c*(b-8*alpha*a)/( (a-b)*(a-c) );
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
 *                = y_{i-1} + 0.5*dx_{i-3/2} * y'_{i-3/2} 
 * where
 *      r_{i-1/2} = y'_{i-1/2} / y'_{i-3/2},
 * and
 *      y'_{i-1/2} = (y_i - y_{i-1})/(x_i - x_{i-1})
 *      y'_{i-3/2} = (y_{i-1} - y_{i-2})/(x_{i-1} - x_{i-2})
 */                        
void AdvectionInterpolationCoefficient::SetFluxLimitedCoefficient(int_t ind, int_t N, const real_t *x, real_t psi, real_t *&deltas){
    real_t x1 = GetXi(x,ind+shiftU1,N);
    real_t x2 = GetXi(x,ind+shiftU2,N);
    real_t dx0 = xf - x1;
    real_t dxf = x1 - x2;
    deltas[2+shiftU1] = 1 + psi * dx0/dxf;
    deltas[2+shiftU2] = -psi*dx0/dxf;
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

/*
1 + .5*a - .5*b
-.5*a
.5*b

a=b=.5:

1
-1/4
.1/4

a=0, b=2
0
0
1

a=2, b=0
2
-1
0
*/

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
                    }
                } else if(bc_lower == AD_BC_DIRICHLET)
                    if(ind==0)
                        for(len_t k=0; k<2*stencil_width; k++)
                            deltas[ir][pind][k] = 0;
 
                if(bc_upper == AD_BC_MIRRORED){
                    for(len_t k=N+stencil_width-ind; k<=k_max; k++){
                        deltas[ir][pind][k_max+2*(N-ind)-k] += deltas[ir][pind][k];
                        deltas[ir][pind][k] = 0;
                    }
                } else if(bc_upper == AD_BC_DIRICHLET)
                    if(ind==N)
                        for(len_t k=0; k<2*stencil_width; k++)
                            deltas[ir][pind][k] = 0;
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
        nnzPerRow = 8*1-1; // = 7
    else
        nnzPerRow = 8*stencil_width-1; // = 15
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
 * Returns the normalized-variable parameter 
 *      PhiHat_{i-1} = (y_{i-1} - y_{i-2}) / (y_i - y_{i-2})
 * which is used to switch between interpolation methods 
 * in order for the flux to satisfy a 'convection boundedness criterion'.
 * See P H Gaskell and A K C Lau, IJNMF 8, 617-641 (1988) for more details
 * on Normalized variable formulations, the convection boundedness criterion
 * and the SMART scheme.
 */
real_t AdvectionInterpolationCoefficient::GetPhiHatNV(int_t ind, int_t N, std::function<real_t(int_t)> y){
    real_t y0 = GetYi(ind+shiftD1, N, y);
    real_t y1 = GetYi(ind+shiftU1, N, y);
    real_t y2 = GetYi(ind+shiftU2, N, y);
   
    if(y0==y2) 
    // returns something smaller than 0 if y1-y2 is negative
    // or something greater than 1 if y1-y2 is positive 
        return 0.5 + 0.6*( (y1>y2) - (y1<=y2) );
    else
        return (y1-y2) / (y0-y2);
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
    return dy0/dy1;
}
