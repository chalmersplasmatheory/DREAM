/**
 * In this file, methods used for the calculation of the bounce averaged
 * delta function appearing in the Rosenbluth-Putvinski avalanche source
 * are collected.
 */





/**
 * Helper function for avalanche deltaHat calculation: 
 * returns the function xi0Star = RESign * sqrt( 1 - 1/BOverBmin * 2/(g+1) )
 * See documentation in doc/notes/theory.
 */
real_t xi0Star(real_t BOverBmin, real_t gamma, int_t RESign){
    real_t p2 = gamma*gamma - 1;
    // This form is numerically stable for arbitrary p and BOverBmin
    return RESign*sqrt( (p2/(gamma+1) + 2*(BOverBmin-1)/BOverBmin ) / (gamma+1) );
}

/**
 * For gsl root finding: returns the function sgn*(xi0Star - xi0), which 
 * defines the integration limits in avalanche deltaHat calculation.
 */
struct xiStarParams {real_t gamma; real_t xi0; len_t ir; real_t Bmin; FluxSurfaceAverager *FSA; int_t sgn; int_t RESign;};
real_t xi0StarRootFunc(real_t theta, void *par){
    struct xiStarParams *params = (struct xiStarParams *) par;

    real_t B, Jacobian, ROverR0, NablaR2;
    params->FSA->GeometricQuantitiesAtTheta(params->ir, theta, B, Jacobian, ROverR0, NablaR2, FLUXGRIDTYPE_DISTRIBUTION);
    real_t Bmin = params->Bmin;
    real_t BOverBmin = 1;
    if(Bmin)
        BOverBmin = B/Bmin;

    return params->sgn*(xi0Star(BOverBmin, params->gamma, params->RESign) - params->xi0);
}


/**
 * Helper function for avalanche deltaHat function: 
 * returns the integrand in the bounce average of the delta function
 * (which has been integrated over a grid cell)
 */
struct hParams {real_t gamma; len_t ir; real_t Bmin; real_t Vp; real_t dxi; FluxSurfaceAverager *FSA; int_t RESign;};
real_t hIntegrand(real_t theta, void *par){
    struct hParams *params = (struct hParams *) par;
    real_t B, Jacobian, ROverR0, NablaR2;
    params->FSA->GeometricQuantitiesAtTheta(params->ir, theta, B, Jacobian, ROverR0, NablaR2, FLUXGRIDTYPE_DISTRIBUTION);

    real_t Bmin = params->Bmin;
    real_t BOverBmin = 1;
    if(Bmin)
        BOverBmin = B/Bmin;

    real_t Vp = params->Vp;
    real_t dxi = params->dxi;
    int_t RESign = params->RESign;

    real_t g = params->gamma;
    real_t xi = RESign*sqrt((g-1)/(g+1));
    real_t xi0 = xi0Star(BOverBmin,g,RESign);
    real_t sqrtgOverP2 = MomentumGrid::evaluatePXiMetricOverP2(xi0,BOverBmin);

    // 2*pi for the trivial phi integral
    return 2*M_PI * xi/xi0 * Jacobian * sqrtgOverP2 / (dxi * Vp);
}


/**
 * This helper function orders the arguments theta1 and theta2  
 * such that an integration from theta1 to theta2 goes via the 
 * high field side (ie doesn't include 0 in the interval)
 */
void orderIntegrationIndicesHFS(real_t *theta1, real_t *theta2){
    if(*theta1 < 0)
        *theta1 += 2*M_PI;
    if(*theta2 < 0)
        *theta2 += 2*M_PI;
    if(*theta2<*theta1){ // ensure that theta1 < theta2
        real_t tmp = *theta1;
        *theta1 = *theta2;
        *theta2 = tmp;
    }

}

/**
 * This helper function orders the arguments theta1 and theta2 
 * such that an integration from theta1 to theta2 goes via the
 * low field side (ie doesn't include pi in the interval)
 */
void orderIntegrationIndicesLFS(real_t *theta1, real_t *theta2){
    if(*theta2<*theta1){ // ensure that theta1 < theta2
        real_t tmp = *theta1;
        *theta1 = *theta2;
        *theta2 = tmp;
    }
    // in order to integrate on low-field side, the smallest theta should be negative
    // if not, the larger theta should be shifted to negative values and be the lower limit
    if(*theta1>0){
        *theta2 -= 2*M_PI;
        real_t tmp = *theta1;
        *theta1 = *theta2;
        *theta2 = tmp;
    }    
}


/**
 * Returns the bounce and cell averaged delta function in xi that
 * appears in the Rosenbluth-Putvinski avalanche source term.
 * 
 * Parameters:
 *       ir: radial grid point
 *        p: momentum
 *     xi_l: xi on the lower cell face
 *     xi_u: xi on the upper cell face
 *       Vp: bounce-integrated metric
 *    VpVol: flux-surface averaged jacobian
 *   RESign: sign of xi of the incident REs (+1 or -1).
 *           Is used to flip the pitch of the source
 */
real_t FluxSurfaceAverager::EvaluateAvalancheDeltaHat(len_t ir, real_t p, real_t xi_l, real_t xi_u, real_t Vp, real_t VpVol, int_t RESign){
    // Since Vp = 0 this point will not contribute to the created density 
    if(Vp==0)
        return 0; //placeholder

    real_t theta_Bmin=0, theta_Bmax=0;
    real_t Bmin = GetBmin(ir, FLUXGRIDTYPE_DISTRIBUTION,&theta_Bmin);
    real_t Bmax = GetBmax(ir, FLUXGRIDTYPE_DISTRIBUTION,&theta_Bmax);
    real_t BmaxOverBmin;

    if(Bmin==Bmax)
        BmaxOverBmin=1;
    else 
        BmaxOverBmin = Bmax/Bmin;
    
    real_t gamma = sqrt(1+p*p);
    if(RESign>=0){
        if( xi0Star(BmaxOverBmin, gamma, RESign) <= xi_l )
            return 0;
        else if( xi0Star(1, gamma, RESign) >= xi_u )
            return 0;
    } else {
        if( xi0Star(1, gamma, RESign) <= xi_l )
            return 0;
        else if( xi0Star(BmaxOverBmin, gamma, RESign) >= xi_u )
            return 0;
    }
    // else, there are two nontrivial intervals [theta_l, theta_u] on which contributions are obtained

    xiStarParams xi_params_u = {gamma,xi_u,ir,Bmin,this, -1, RESign}; 
    xiStarParams xi_params_l = {gamma,xi_l,ir,Bmin,this, 1, RESign}; 

    // from now on the logic in this function is a proper mind fuck, 
    // and I apologize to future maintainers who have to touch this. 
    // Godspeed.
    //                                 -- Ola, Uddevalla, 2020-10-07

    // check if xi0Star < xi_u is satisfied for all theta
    bool upperForAllTheta = (xi0StarRootFunc(theta_Bmin, &xi_params_u) > 0) && (xi0StarRootFunc(theta_Bmax, &xi_params_u) > 0);
    // check if xi0Star > xi_l is satisfied for all theta
    bool lowerForAllTheta = (xi0StarRootFunc(theta_Bmin, &xi_params_l) > 0) && (xi0StarRootFunc(theta_Bmax, &xi_params_l) > 0);

    // if all poloidal angles contribute fully to the integral, return the known exact value.
    if(upperForAllTheta && lowerForAllTheta)
        return 2*M_PI*VpVol/(Vp/(p*p)) *  rGrid->GetFSA_B(ir) / (p*p*(xi_u-xi_l));


    hParams h_params = {gamma,ir,Bmin,Vp,xi_u-xi_l, this, RESign};
    gsl_function h_gsl_func;
    h_gsl_func.function = &(hIntegrand);
    h_gsl_func.params = &h_params;
    
    // settings for integral
    real_t epsabs = 0, epsrel = 1e-4, lim = gsl_adaptive->limit, error;
    real_t deltaHat;


    gsl_function gsl_func;
    gsl_func.function = &(xi0StarRootFunc);
    real_t theta_u1, theta_u2, theta_l1, theta_l2;


    // if below is satisfied, integrate between theta_l1 and theta_l2 via 
    // the poloidal angles where the inequalities (in gsl_func) are satisfied
    if(upperForAllTheta){
        gsl_func.params = &xi_params_l;
        FindThetas(theta_Bmin,theta_Bmax,&theta_l1, &theta_l2, gsl_func, gsl_fsolver);

        if(RESign==1)
            orderIntegrationIndicesHFS(&theta_l1,&theta_l2);
        if(RESign==-1)
            orderIntegrationIndicesLFS(&theta_l1,&theta_l2);
        
        real_t pts[2] = {theta_l1,theta_l2};
        int npts = 2;
        gsl_integration_qagp(&h_gsl_func,pts,npts,epsabs,epsrel,lim,gsl_adaptive,&deltaHat, &error);
        return deltaHat;
    }

    // like previous block for theta_u1 and theta_u2
    if(lowerForAllTheta){
        gsl_func.params = &xi_params_u;
        FindThetas(theta_Bmin,theta_Bmax,&theta_u1, &theta_u2, gsl_func, gsl_fsolver);
        if(RESign==1)
            orderIntegrationIndicesLFS(&theta_u1,&theta_u2);
        if(RESign==-1)
            orderIntegrationIndicesHFS(&theta_u1,&theta_u2);

        real_t pts[2] = {theta_u1,theta_u2};
        int npts = 2;
        gsl_integration_qagp(&h_gsl_func,pts,npts,epsabs,epsrel,lim,gsl_adaptive,&deltaHat, &error);
        return deltaHat;        
    }

    // otherwise, integrate between theta_u1 and theta_l1 and between theta_u2 and theta_l2.
    // These should be ordered such that the intervals do not cross theta_Bmax or theta_Bmin
    gsl_func.params = &xi_params_u;
    FindThetas(theta_Bmin,theta_Bmax,&theta_u1, &theta_u2, gsl_func, gsl_fsolver);
    
    gsl_func.params = &xi_params_l;
    FindThetas(theta_Bmin,theta_Bmax,&theta_l1, &theta_l2, gsl_func, gsl_fsolver);

    real_t deltaHat1, deltaHat2;
    int npts = 2;
    real_t pts1[2] = { min(theta_l1,theta_u1), max(theta_l1,theta_u1) };
    real_t pts2[2] = { min(theta_l2,theta_u2), max(theta_l2,theta_u2) };
    gsl_integration_qagp(&h_gsl_func,pts1,npts,epsabs,epsrel,lim,gsl_adaptive,&deltaHat1, &error);
    gsl_integration_qagp(&h_gsl_func,pts2,npts,epsabs,epsrel,lim,gsl_adaptive,&deltaHat2, &error);

    return deltaHat1+deltaHat2;
}
