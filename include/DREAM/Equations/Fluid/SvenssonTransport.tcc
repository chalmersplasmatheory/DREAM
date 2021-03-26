/**
 * Implementation of a diffusive transport term which can be applied to both
 * kinetic and fluid grids, and which allows one to prescribe the diffusion
 * coefficient in time and phase space.
 */

#include <type_traits>
#include "DREAM/Equations/Fluid/SvenssonTransport.hpp"

/**
 * Constructor.
 */
template<typename T>
DREAM::SvenssonTransport<T>::SvenssonTransport(
    DREAM::FVM::Grid *grid,
    real_t pStar,
    DREAM::FVM::UnknownQuantityHandler *unknowns,
    DREAM::RunawayFluid *REFluid,
    struct dream_4d_data *inputData4dStruct,
    enum FVM::Interpolator1D::interp_method timeInterpMethod 
) : T(grid),
    nr_f(grid->GetNr()+1),
    nt(inputData4dStruct->nt), nr(inputData4dStruct->nr),
    np1(inputData4dStruct->np1), np2(inputData4dStruct->np2),
    np(CountNp(np1,pStar,inputData4dStruct->p1)), nxi(CountNxi(np2)),
    EID(unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD)),
    pStar(pStar),
    t(inputData4dStruct->t), r(inputData4dStruct->r),
    p1(inputData4dStruct->p1), p2(inputData4dStruct->p2),
    coeff4dInput(inputData4dStruct->x), // Size nt-by-(nr*np2*np1)
    inputMomentumGridType(inputData4dStruct->gridtype),
    inputInterp3dMethod(inputData4dStruct->ps_interp),
    unknowns(unknowns), REFluid(REFluid),
    timeInterpMethod(timeInterpMethod)
{
    // Checks that pStar is valid.
    if ( pStar <= 0 ) {
        throw DREAMException(
            "Invalid lower momentum bound, pStar=%.2f, must be >0.",
            pStar
            );
    }
    // Checks that the input data is defined on a p-xi momentum grid.
    if (inputMomentumGridType != FVM::Interpolator3D::GRID_PXI ){
        // SvenssonTransport is, at the moment, only defined for p-xi.
        throw DREAMException(
            "Wrong input momentum grid, input must be specified on p-xi grid."
            );
    }
    else{
        xi=p2;
        SetMomentumCoordinate();
    }


    this->coeffTRXiP = new real_t[ nt * nr_f * np * nxi ];
    this->coeffRP = new real_t[ nr_f * np ];
    this->integrand = new real_t[ np ];

    InterpolateCoefficient();
}


/**
 * Destructor.
 */
template<typename T>
DREAM::SvenssonTransport<T>::~SvenssonTransport() {
    delete [] this->t;
    delete [] this->r;
    delete [] this->p1;
    delete [] this->p2;
    delete [] this->coeff4dInput;
    //delete [] this->xi;// Should xi=p2 be deleted the same way as p2
    delete [] this->p;
    delete [] this->coeffTRXiP;
    delete [] this->coeffRP;
    delete [] this->integrand;
    if (this->interpTCoeff != nullptr) {
        delete this->interpTCoeff;
    }
}


/**
 * Interpolate the time-depentent input coefficient onto the DREAM
 * radial flux-grid, then storing that data in `interpTCoeff`, which
 * is then used to evaluate the time dependance of the coefficients.
 *
 * This function creates a 1D interpolator, which is then used to
 * evaluate the time dependence of the coefficients.
 */
template<typename T>
void DREAM::SvenssonTransport<T>::InterpolateCoefficient() {    
    const len_t N  =  this->nr_f * this->nxi * this->np;
    
    for (len_t it = 0, offset = 0; it < nt; it++) {
        DREAM::FVM::Interpolator3D intp3d_tmp(
            nr, np2, np1, r, p2, p1, coeff4dInput[it],
            inputMomentumGridType, inputInterp3dMethod, false
            );
        // Interpolating the coefficients (of every supplied time
        // step) onto the r_f grid used by DREAM.
        intp3d_tmp.Eval(nr_f, nxi, np, this->grid->GetRadialGrid()->GetR_f(), xi,p,
                        FVM::Interpolator3D::momentumgrid_type::GRID_PXI, coeffTRXiP+offset);
        // This is more memory intese, than doing the interpolation
        // onto the r_f grid in the xiAverage function. However, this
        // method gives faster simulation runtimes due to otherwise
        // having to redo many of the r_f interpolations at every time
        // step.
        offset+=N;
    }
    // Note that `coeffTRXiP` now contains r_f, xi and p data for _every_ time step.

    if (this->interpTCoeff != nullptr) {
        delete this->interpTCoeff;
    }
    this->interpTCoeff = new DREAM::FVM::Interpolator1D(nt, N, t, coeffTRXiP, timeInterpMethod);
}




/**
 * Rebuild this term by evaluating and setting the advection and
 * diffusion terms for the next time step.
 */
template<typename T>
void DREAM::SvenssonTransport<T>::Rebuild(
    const real_t time, const real_t, DREAM::FVM::UnknownQuantityHandler* 
    ) {
    // Note that the Interp1D object doesn't allocate new memory in
    // this process (with "nearest" interpolation method).
    const real_t *coeffRXiP = this->interpTCoeff->Eval(time);

    
    // // Printing out the radially interpolated coefficients
    // for (len_t ir=0; ir<nr_f; ir++){
    //     printf("ir = %3lu | r = %0.3f\n", ir, this->grid->GetRadialGrid()->GetR_f(ir));
    //     for(len_t i=0; i<np; i++){
    //         printf("i  = %2lu  | p = %03.3f   ", i,p[i]);
    //         for(len_t j=0; j<nxi; j++){
    //             printf("%+3.3f ",coeffRXiP[(ir*nxi+j)*np+i]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");fflush(stdout);
    // }
    

    // Note that the xi-average has to be redone for every timestep,
    // since the prescribed distribution function changes with time.
    xiAverage(coeffRXiP);

    // Make sure that there are enough p-points to do a trapz integration.
    if (np > 1){
        #define COEFFDEBUG false // DEBUG flag
        real_t coeff_abs_max = 0.0; // DEBUG
        len_t i_max=0;
        // Iterate over the radial flux grid
        for (len_t ir = 0; ir < this->nr_f; ir++) {
            
            // The varaible to be added to
            real_t pIntCoeff = 0;
            
            this->EvaluateIntegrand(ir);
            
            // The trapz integration in p
            for (len_t i = 0; i < this->np-1; i++) {
                pIntCoeff += 0.5 * (this->integrand[i] + this->integrand[i+1])
                    * (this->p[i+1] - this->p[i]);
                // YYY Jacobian??? * this->grid->GetVp(ir,i,0); 
            }
            this->_setcoeff(ir, pIntCoeff);
            if(COEFFDEBUG){
                if ( abs(pIntCoeff) > abs(coeff_abs_max) ){
                    coeff_abs_max = pIntCoeff;
                    i_max=ir;
                }
            }
        }
        if(COEFFDEBUG){
            printf("%s", typeid(T).name());
            printf(" | coef_abs_max = %0.3e @ ir = %3lu\n", coeff_abs_max,i_max);
            fflush(stdout); // DEBUG
        }
    }
    else{
        for (len_t ir = 0; ir < this->nr_f; ir++) {
            // Sets the coefficient in the event that the integration
            // cannot be performed. (We should never get here though!)
            //
            // `coeffRP` is here of size nr_f*np = nr_f*1
            this->_setcoeff(ir, coeffRP[ir]);
        }
    }
}






/**
 * Function for evaluating the xi average of the coefficient over the
 * prescribed pitch angle distribution function.
 */
template<typename T>
void DREAM::SvenssonTransport<T>::xiAverage(const real_t *coeffRXiP){
    // Input `coeffRXiP` is the r-xi-p coefficient data (nr_f * nxi * np)
    // Writing xi-averaged data to `this->coeffRP` (nr_f * np)

    #define XIDEBUG false // DEBUG flag

    if (XIDEBUG){
        printf("         |");                                // DEBUG
        printf("           |");                              // DEBUG
        for(len_t i=0; i<np;i++) printf("  ip |  coeff |");  // DEBUG
        printf("\n");                                        // DEBUG
    }

    // Make sure that there are enough xi-points to do a trapz integration.
    if (nxi > 1) {
        // Arrays with the E field, collision times and Zeff radial profiles
        const real_t *EVec    = this->unknowns->GetUnknownData(this->EID);
        const real_t *tauVec  = this->REFluid->GetElectronCollisionTimeRelativistic();
        const real_t *ZeffVec = this->REFluid->GetIonHandler()->GetZeff(); 

        // Normalization factor  fot the E field
        const real_t normFactor = Constants::ec / (Constants::me * Constants::c );
        
        for (len_t ir=0, offset=0; ir < nr_f ; ir++){ 
            if(XIDEBUG){
                printf("ir = %3lu | ",ir); // DEBUG
                printf("r = %0.3f | ",this->grid->GetRadialGrid()->GetR_f()[ir]); // DEBUG
            }

            // Evaluate E and Zeff on the fluxgrid.
            real_t Zeff_f = this->EvalOnFluxGrid(ir,ZeffVec); // local Z_eff
            real_t tau_f = this->EvalOnFluxGrid(ir, tauVec); // local collision time
            real_t E_f = this->EvalOnFluxGrid(ir, EVec);    // local E field
                                   

            for (len_t i=0; i < this->np ; i++){ 
                
                
                real_t p_sq = this->p[i] * this->p[i]; // p^2
                // Calculating the pitch-angle width, w, of the
                // prescribed distribution.
                real_t w = 2.0 * normFactor * abs(E_f) * p_sq * tau_f
                    / ( (1.+Zeff_f) * sqrt(1+p_sq) );
                // YYY Think about if abs(E_f) affects advection
                // Check with Konsta

                // Helper variable for the integrand prefactor.
                // Factor 0.5 comes from trapz integration method
                real_t w_factor = 0.5 * w / (1.0 - exp( -2.0 * w ) );

                // The relevant indices are: (ir*nxi+j)*np+i and the same with j+1.
                len_t ind = offset*nxi + i;
                // Variables containing the integrands
                real_t g1 = coeffRXiP[ind] * exp( w*(this->xi[0]-1.0) );// init first integrand
                real_t g2;
                

                real_t avg = 0; // Varaible containing the xi average

                //Trapz integration in xi.
                for (len_t j=0; j<nxi-1; j++){
                    // YYY Possibly change this trapz method to a GSL
                    // integral.
        
                    // The relevant indices are: (ir*nxi+j)*np+i and the same with j+1.
                    // g1 is calculated in previous interation.
                    // We want the next index for g2:
                    ind+=np;
                    g2 = coeffRXiP[ind] * exp( w * (this->xi[j+1]-1.0) );

                    
                    // N.B! a Prefactor 0.5 has been included in w_factor
                    avg += w_factor * (g1 + g2) * (xi[j+1]-xi[j]);

                    // The next iteration's g1 is the current one's g2:
                    g1=g2;
                }
                this->coeffRP[i+offset] = avg / (xi[nxi-1]-xi[0]);
                if(XIDEBUG)
                    printf("%3lu | %0.4f | ", i, this->coeffRP[i+offset]); // DEBUG
            }
            offset+=this->np;
            if(XIDEBUG) printf("\n"); fflush(stdout); // DEBUG
        }
    }
    else{
        // If there is only one xi-point, we simply take the
        // coefficient value at face value, no modifications.
        for (len_t ir=0, offset=0; ir < nr_f ; ir++){
            for (len_t i=0; i < this->np ; i++){
                // `coeffRXiP` is here of size nr_f*nxi*np = nr_f*1*np
                this->coeffRP[i+offset] = coeffRXiP[offset + i];
            }
            offset+=this->np;
        }
    }
    if(XIDEBUG) printf("Xi-averaging, done!\n"); fflush(stdout); // DEBUG
    #undef XIDEBUG
}




/**
 * Helper function for calculating the inverse of p-bar, with the
 * (optional) additional calculation of the derivative of
 * p-bar-inverse. 
 *
 * `p-bar` is the name given to the factor dividing `-(p - p*)` in the
 * exponential of eqn (4.2) in Svensson et al. 2020
 * [https://arxiv.org/abs/2010.07156v1].
 * 
 * These values are calculated on the flux grid, meaning that
 * interpolation (and extrapolation) from the cell grid is being
 * performed. This is done via inter-/extrapolation of p-bar-inverse,
 * instead of first inter-/extrsapolating the values going into p-bar.
 *
 * YYY XXX Consider using `EvalOnFluxGrid` for much of the
 * functionality of this function instaead. However, derivative will
 * be tricky to replace!
 */
template<typename T>
real_t DREAM::SvenssonTransport<T>::GetPBarInv_f(len_t ir, real_t *dr_pBarInv_f){
    // Need interpolation from cell grid to flux grid:
    // pBar_f[0]           = pBar[0]
    // pBar_f[0<ir<nr_f-1] = interpolate
    // pBar_f[nr_f-1]      = extrapolate

    // Inverse of p-bar on the Flux grid, with additional helper variable.
    real_t pBarInv_f, tmp_pBarInv_f; 

    // Essential values taken on the raidal (cell) grid
    const real_t *E       = this->unknowns->GetUnknownData(this->EID);
    const real_t *EcEff   = this->REFluid->GetEffectiveCriticalField();
    const real_t *gamma_r = this->REFluid->GetAvalancheGrowthRate();

    // Normalization factor for pBarInv
    const real_t normFactor = Constants::me * Constants::c / Constants::ec;

    // Grid step size in the radial grid for the derivative.
    const real_t *dr_f    = this->grid->GetRadialGrid()->GetDr_f();
    const real_t *dr      = this->grid->GetRadialGrid()->GetDr(); 
    
    // Interpolating (extrapolating) the inverse of p-bar onto the
    // flux grid.
    if (ir == 0) {
        // Zero flux at r = 0. Therefore choose the value at "ir=1/2".
        pBarInv_f = normFactor * gamma_r[0] / (std::abs(E[0])-EcEff[0]);
        
        if(dr_pBarInv_f != nullptr)
            *dr_pBarInv_f = 0.0;
    }
    else if (ir == this->nr_f - 1) {
        // Linearly extrapolating the value at the end point from the
        // two previous points.
        pBarInv_f     = normFactor * gamma_r[ir-1] / (std::abs(E[ir-1])-EcEff[ir-1]);
        tmp_pBarInv_f = normFactor * gamma_r[ir-2] / (std::abs(E[ir-2])-EcEff[ir-2]);

        if(dr_pBarInv_f != nullptr){
            // N.B.! This order of operations is important
            // Derivative:
            *dr_pBarInv_f = (pBarInv_f - tmp_pBarInv_f) / dr_f[ir-2];
            // Extrapolation:
            pBarInv_f += (*dr_pBarInv_f) * 0.5*dr[ir-1];
        }
        else{
        pBarInv_f += (pBarInv_f - tmp_pBarInv_f) * 0.5*dr[ir-1]/dr_f[ir-2];
        }
        // The above is the same as:
        // pBarInv_f = tmp_pBarInv_f
        //     + (pBarInv_f - tmp_pBarInv_f)
        //     * (0.5*dr[ir-2] + dr[ir-1])/dr_f[ir-2];
    }
    else {
        // In the middle, we simply linearly interpolate
        tmp_pBarInv_f = normFactor * gamma_r[ir-1] / (std::abs(E[ir-1])-EcEff[ir-1]);
        pBarInv_f  = normFactor * gamma_r[ir] / (std::abs(E[ir])-EcEff[ir]);

        // N.B.! This order of operations is important!
        if(dr_pBarInv_f != nullptr){
            // Derivative:
            *dr_pBarInv_f = (pBarInv_f - tmp_pBarInv_f) / dr_f[ir-1]; 
            // Interpolation:
            pBarInv_f = tmp_pBarInv_f + (*dr_pBarInv_f) * 0.5*dr[ir-1];
            //pBarInv_f -= (*dr_pBarInv_f) * 0.5*dr[ir];
        }
        else{
            //pBarInv_f = ( dr[ir]*tmp_pBarInv_f + dr[ir-1]*pBarInv_f ) * 0.5 / dr_f[ir-1];
            pBarInv_f = tmp_pBarInv_f + (pBarInv_f - tmp_pBarInv_f) * 0.5*dr[ir-1]/dr_f[ir-1];
        }
        // This is for uniform step sizes:
        // pBarInv_f += tmp_pBarInv_f;
        // pBarInv_f *= 0.5;
    }

    return pBarInv_f;
}




/**
 * Function for evalutaing a vector, defined in the cell grid, on the
 * flux grid. Interpolating in the inner, extrapolating on the outer
 * end, and setting zero flux on the inner most point.
 */ 
template<typename T>
real_t DREAM::SvenssonTransport<T>::EvalOnFluxGrid(len_t ir, const real_t *vec){
    // Need interpolation from cell grid to flux grid:
    // vec_f[0]         = vec[0]
    // vec_f[0<ir<nr_f] = interpolate
    // vec_f[nr_f]      = extrapolate


    // The value of the vec quantity on the flux grid
    real_t vec_f;

    // Grid step size in the radial grid for the derivative.
    const real_t *dr_f = this->grid->GetRadialGrid()->GetDr_f();
    const real_t *dr = this->grid->GetRadialGrid()->GetDr(); 

    
    // Interpolating (extrapolating) the vec onto the flux grid.
    if (ir == 0) {
        // Zero flux at r = 0. Therefore choose the value at "ir=1/2".
        vec_f = vec[0];
    }
    else if (ir == this->nr_f - 1) {
        // Linearly extrapolating the value at the end point from the
        // two previous points.
        vec_f = vec[ir-1] + (vec[ir-1] - vec[ir-2]) * 0.5*dr[ir-1]/dr_f[ir-2];
    }
    else {
        // In the middle, we simply linearly interpolate
        vec_f = vec[ir-1] + (vec[ir] - vec[ir-1]) * 0.5*dr[ir-1]/dr_f[ir-1];
    }

    return vec_f;
}



/**
 * Counting the number of useable p points.
 */
template<typename T>
const len_t DREAM::SvenssonTransport<T>::CountNp(
    const len_t np1In, const real_t pStarIn, const real_t *p1In
    ) {
    len_t np_count=0;
    for (len_t i1 = 0; i1 < np1In; i1++){
        // Check that p1 is strictly increasing.
        if (i1 > 0){
            if (p1In[i1] <= p1In[i1-1])
                throw DREAMException("Input momentum coordiantes must be strictly increasing.");
        }
        // Add one to the count for every value larger than pStar.
        np_count += (p1In[i1]>pStarIn ? 1 : 0 );
    }

    // If pStar is within the range of p1, we also include pStar in p.
    np_count += (p1In[0]>pStarIn ? 0 : 1 );

    // There must be at least two valid p values
    if ( np_count  < 2 )
        throw DREAMException("Not enough coordiantes p>pStar supplied.");
    
    return np_count;
}

/**
 * Setting the p array with only the valid p coordinates
 * (p1[i]>pStar).
 */
template<typename T>
void DREAM::SvenssonTransport<T>::SetMomentumCoordinate() {
    p = new real_t[np];
    // Offset to the initilization index if we include pStar or not,
    // wich depends on whether or not pStar is within the range of p1.
    len_t pStarOffset = (p1[0]>pStar ? 0 : 1 );
    // Set p[0] to pStar. If pStarOffset==0 (i.e. don't include
    // pStar), then p[0] will be overwritten below.
    p[0] = pStar;

    // Set the rest of the valid p values.
    for (len_t ip = pStarOffset; ip < np; ip++){
        p[ip]=p1[np1+ip-np];
    }
}
