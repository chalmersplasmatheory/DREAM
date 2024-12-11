/**
 * In this file, methods used for the calculation of the bounce averaged
 * Chiu-Harvey avalanche source are collected.
 */

/**
 * Helper function for avalanche bounce average calculation: 
 * returns the function 
 * xi0_min = RESign * sqrt((gamma - 1) / (gamma + 1) 
 *                          * (gamma_{in,max} + 1) / (gamma_{in,max} - 1))
 */
real_t ximin(real_t gamma, real_t gamma_max/*, int_t RESign*/){
    return /*RESign**/sqrt((gamma - 1) / (gamma + 1) * (gamma_max + 1) / (gamma_max - 1));
}

/**
 * Helper function for avalanche bounce average calculation: 
 * returns the function 
 * xi0_max = RESign * sqrt(gamma / (gamma + 1))
 */
real_t ximax(real_t gamma/*, int_t RESign*/){
    return /*RESign**/sqrt(gamma / (gamma + 1));
}

/*
 * Integrand of each bounce average calculation term.
 */
real_t FluxSurfaceAverager::BA_CH(real_t xiOverXi0, real_t BOverBmin, real_t, real_t, void *par){
    struct avParams *params = (struct avParams *) par;
    real_t xi0 = params->xi0;
    real_t gamma = params->gamma;
    real_t gamma_max = params->gamma_max;
    int_t iTerm = params->iTerm;
    real_t BminOverBmax = params->BminOverBmax;
    
    if (ximax(gamma/*, 1*/) < ximin(gamma, gamma_max/*, 1*/))
        return 0;
    if (ximax(gamma/*, 1*/) < xiOverXi0*abs(xi0))
        return 0;
    if (ximin(gamma, gamma_max/*, 1*/) > xiOverXi0*abs(xi0))
        return 0;

    real_t xi2 = 1. - (1. - xi0*xi0) * BOverBmin;
    real_t denom = xi2 * (gamma + 1.) - (gamma - 1.);
    real_t term; 
    switch(iTerm){
        case 0:
            term = xiOverXi0 / (BOverBmin*BOverBmin) / (denom*denom*denom*denom);
            break; 
        case 1:
            term = xiOverXi0 * xi2 / (BOverBmin*BOverBmin) / (denom*denom*denom*denom);
            break; 
        case 2:
            term = xiOverXi0 * xi2*xi2 / (BOverBmin*BOverBmin) / (denom*denom*denom*denom);
            break; 
        case 3:
            term = xiOverXi0 / BOverBmin / (denom*denom*denom);
            break; 
        case 4:
            term = xiOverXi0 * xi2 / BOverBmin / (denom*denom*denom);
            break; 
        case 5:
            term = xiOverXi0 * xi2*xi2 / BOverBmin / (denom*denom*denom);
            break; 
        default:
            throw FVMException("There is no '%d'th term in bounce average of Chiu-Harvey avalanche operator.", term);
    }
    if ((1-xi0*xi0) > BminOverBmax)
        return 0.5*term;
    return term;
}

/*
 * Integrand from discretization in xi.
 */
real_t FluxSurfaceAverager::integrandXi(real_t xi0, void *par){
    struct intXiParams *params = (struct intXiParams *) par;
    len_t ir = params->ir;
    fluxGridType fgt = params->fgt;
    real_t gamma = params->gamma;
    real_t gamma_max = params->gamma_max;
    real_t BminOverBmax = params->BminOverBmax;
    //int_t RESign = params->RESign;
    int_t iTerm = params->iTerm;
    FluxSurfaceAverager *FSA = params->FSA;
    if (xi0 < 0.)
        return 0;
    if(xi0 < ximin(gamma, gamma_max/*, RESign*/))
        return 0;
    /*if ((RESign>=0 && xi0 < 0.) || (RESign<0 && xi0 > 0.)) { 
        return 0.;
    }
    if(RESign>=0){
        if(xi0 < ximin(gamma, gamma_max, RESign))
            return 0;
    } else {
        if(xi0 > ximin(gamma, gamma_max, RESign))
            return 0;
    }*/
    
    avParams avg_params[5] = {xi0, gamma, gamma_max, iTerm, BminOverBmax};
    real_t factor_BA = FSA->CalculatePXiBounceAverageAtP(ir, xi0, fgt, &BA_CH, avg_params);
    
    if (iTerm < 3)
        return xi0*xi0 / ((1 - xi0*xi0) * (1 - xi0*xi0)) * factor_BA;
    return xi0*xi0 / (1 - xi0*xi0) * factor_BA;
}

/*
 * Integrand from discretization in momentum.
 */
real_t FluxSurfaceAverager::integrandP(real_t p, real_t gamma_max, len_t ir, real_t xi_l, real_t xi_u, real_t BminOverBmax, fluxGridType fgt/*, int_t RESign*/){
    real_t gamma = sqrt(1+p*p);
    real_t epsabs = 0, epsrel = 1e-8, lim = gsl_ws_CH->limit, error;
    
    gsl_function int_gsl_func;
    int_gsl_func.function = &(integrandXi);

    real_t terms[6]; 
    for(int_t i=1; i<6; i++){
        intXiParams intXi_params = {ir, fgt, gamma, gamma_max, BminOverBmax, /*RESign, */i, this};
        int_gsl_func.params = &intXi_params;
        gsl_integration_qag(&int_gsl_func,xi_l,xi_u,epsabs,epsrel,lim,QAG_KEY,gsl_ws_CH,&terms[i], &error);
    }



    return 1/(gamma*p)*
            (4*(gamma-1)/(gamma+1)*terms[0]
              +8*terms[1]
              +4*(gamma+1)/(gamma-1)*terms[2]
              -1*(gamma-1)*(gamma*gamma-2)*terms[3]
              +2*(gamma+1)*(gamma*gamma-gamma-3)*terms[4]
              -(gamma+1)*(gamma+1)*(gamma*gamma-2*gamma+4)/(gamma-1)*terms[5]);
}
/*
// TODO: Remove
real_t FluxSurfaceAverager::FSA_CH(real_t BOverBmin, real_t, real_t, void *par){
    struct avParams *params = (struct avParams *) par;
    real_t xi0 = params->xi0;
    real_t gamma = params->gamma;
    int_t iTerm = params->iTerm;
    int_t RESign = params->RESign;
    
    if (RESign >= 0) {
        if (xiOverXi0 <= 0)
            return 0;
        if (ximax(gamma, RESign) < xiOverXi0*abs(xi0))
            return 0;
    } else {
        if (xiOverXi0 >= 0)
            return 0;
        if (ximax(gamma, RESign) > xiOverXi0*abs(xi0))
            return 0;
    }

    real_t xi2 = 1. - (1. - xi0*xi0) * BOverBmin;
    real_t denom = xi2 * (gamma + 1.) - (gamma - 1.);
    real_t term; 
    switch(iTerm){
        case 0:
            term = 1 / BOverBmin / (denom*denom*denom*denom);
        case 1:
            term = xi2 / BOverBmin / (denom*denom*denom*denom);
        case 2:
            term = xi2*xi2 / BOverBmin / (denom*denom*denom*denom);
        case 3:
            term = xiOverXi0 / (denom*denom*denom);
        case 4:
            term = xi2 / (denom*denom*denom);
        case 5:
            term = xi2*xi2 / (denom*denom*denom);
        default:
            throw FVMException("There is no '%d'th term in bounce average of Chiu-Harvey avalanche operator.", term);
    }
    return term;
}

real_t FluxSurfaceAverager::integrandXi_passing_test(real_t xi0, void *par){
    if (xi0 < 0.) {
        return 0.;
    }
    struct intXiParams *params = (struct intXiParams *) par;
    len_t ir = params->ir;
    fluxGridType fgt = params->fgt;
    real_t gamma = params->gamma;
    //real_t gamma_max = params->gamma_max;
    int_t RESign = params->RESign;
    int_t iTerm = params->iTerm;
    FluxSurfaceAverager *FSA = params->FSA;

    if(RESign>=0){
        if(xi0 < ximin(gamma, RESign))
            return 0;
    } else {
        if(xi0 > ximin(gamma, RESign))
            return 0;
    }

    avParams avg_params = {xi0, gamma, iTerm};
    real_t factor_FSA = FSA->SurfaceAverage(ir, fgt, &FSA_CH, avg_params);
    
    if (iTerm < 3)
        return xi0 / ((1 - xi0*xi0) * (1 - xi0*xi0)) * factor_BA;
    return xi0 / (1 - xi0*xi0) * factor_BA;
}

real_t FluxSurfaceAverager::integrandP_passing_test(real_t p, void *par){
    struct intPParams *params = (struct intPParams *) par;
    len_t ir = params->ir;
    real_t xi_l = params->xi_l;
    real_t xi_u = params->xi_u;
    fluxGridType fgt = params->fgt;
    //real_t p_max = params->gamma_max;
    int_t RESign = params->RESign;
    FluxSurfaceAverager *FSA = params->FSA;

    real_t gamma = sqrt(1+p*p);
    //real_t gamma_max = sqrt(1+p_max*p_max);
    real_t epsabs = 0, epsrel = 1e-4, lim = gsl_ws_CH->limit, error;
    
    gsl_function int_gsl_func;
    int_gsl_func.function = &(integrandXi_passing_test);

    real_t terms[6]; 
    for(int_t i=1; i<6; i++){
        intXiParams intXi_params = {ir, fgt, gamma, RESign, i, FSA};
        gsl_integration_qag(&int_gsl_func,xi_l,xi_u,epsabs,epsrel,lim,QAG_KEY,gsl_ws_CH,&terms[i], &error);
    }

    return p/gamma*
            (4*(gamma-1)/(gamma+1)*terms[0]
              +8*terms[1]
              +4*(gamma+1)/(gamma-1)*terms[2]
              -1*(gamma-1)*(gamma*gamma-2)*terms[3]
              +2*(gamma+1)*(gamma*gamma-gamma-3)*terms[4]
              -(gamma+1)*(gamma+1)*(gamma*gamma-2*gamma+4)/(gamma-1)*terms[5]);
}
*/

// TODO: Ta bort Vp, VpVol

/**
 * Returns the bounce and cell averaged Chiu-Harvey avalanche
 * source term.
 * 
 * Parameters:
 *       ir: radial grid point
 *      p_l: momentum on the lower cell face
 *      p_u: momentum on the upper cell face
 *    p_max: maxmimum momentum of RE grid TODO: Remove
 *     xi_l: xi on the lower cell face
 *     xi_u: xi on the upper cell face
 *   RESign: sign of xi of the incident REs (+1 or -1).
 *           Is used to flip the pitch of the source
 *      fgt: fluxGridType object needed for bounce averaging.
 */
real_t FluxSurfaceAverager::EvaluateAvalancheCHBounceAverage(len_t ir, real_t p_i, real_t p_max, real_t xi_l, real_t xi_u,  fluxGridType fgt/*real_t Vp, real_t VpVol,*/){
    real_t theta_Bmin=0, theta_Bmax=0;
    real_t Bmin = GetBmin(ir, FLUXGRIDTYPE_DISTRIBUTION,&theta_Bmin);
    real_t Bmax = GetBmax(ir, FLUXGRIDTYPE_DISTRIBUTION,&theta_Bmax);
    real_t BminOverBmax;
    real_t xi_u_max, xi_l_max; 
    real_t Delta_xi = xi_u - xi_l;
    real_t xi_j = (xi_u + xi_l) / 2.;
    real_t gamma_i = sqrt(1+p_i*p_i);
    real_t gamma_max = sqrt(1+p_max*p_max);
    //printf("\nximin=%.4f, ximax=%.4f", ximin(gamma_i, gamma_max/*, 1*/), ximax(gamma_i/*, 1*/));
    //printf("\np=%.4f, gamma=%.4f, xi=%.4f, xi_l=%.4f, xi_u=%.4f, p_max=%.4f, gamma_max=%.4f\n",p_i, gamma_i, xi_j, xi_l, xi_u, p_max, gamma_max);
    
    // TODO: Ok?
    //int_t RESign = 1;
    if (xi_j < 0){
        real_t xi_l_temp = xi_l;
        xi_l = -xi_u;
        xi_u = -xi_l_temp;
        xi_j = -xi_j;
        //RESign = -1;
    }
    // TODO: Dubbetänk här: var ska det vara max och var ska det vara min (0 (ingen label))?
    if(Bmin==Bmax) {
        BminOverBmax=1;
        xi_u_max = xi_u; 
        xi_l_max = xi_l;
    } else {
        BminOverBmax = Bmin/Bmax;
        xi_u_max = sqrt(1 - (1 - xi_u*xi_u)/BminOverBmax);
        if (xi_u < 0)
            xi_u_max *= -1.;
        xi_l_max = sqrt(1 - (1 - xi_l*xi_l)/BminOverBmax);
        if (xi_l < 0)
            xi_l_max *= -1.;
    }
    //if(RESign>=0){
        if( ximax(gamma_i/*, RESign*/) <= xi_l )
            return 0;
        else if( ximin(gamma_i, gamma_max/*, RESign*/) >= xi_u_max )
            return 0;
    /*} else {
        if( ximin(gamma_i, gamma_max, RESign) <= xi_l_max )
            return 0;
        else if( ximax(gamma_i, RESign) >= xi_u )
            return 0;
    }*/
    
    /*
    real_t epsabs = 0, epsrel = 1e-4, lim = gsl_ws_CH->limit, error;

    gsl_function int_gsl_func;
    int_gsl_func.function = &(integrandP);

    real_t avaCH_BA; 
    intPParams intP_params = {ir, xi_l, xi_u, fgt, RESign, this, gsl_ws_CH, QAG_KEY};
    gsl_integration_qag(&int_gsl_func,p_l,p_u,epsabs,epsrel,lim,QAG_KEY,gsl_ws_CH,&avaCH_BA, &error);
    */
    real_t avaCH_BA = 1 / (abs(xi_j) * Delta_xi) * integrandP(p_i, gamma_max, ir, xi_l, xi_u, BminOverBmax, fgt/*, RESign*/);
    printf("BA=%.8e\n", avaCH_BA);
    if (avaCH_BA < 0)// TODO: remove
        printf("\nNegative BA! BA=%.4e, xi_j=%.4e\n", avaCH_BA, xi_j);
    /*
    if  ( (1-xi0*xi0) >  BminOverBmax){
        gsl_function int_gsl_func;
        int_gsl_func.function = &(integrandP_passing_test);

        real_t avaCH_FVA_passing_test; 
        intPParams intP_params = {ir, fgt, xi_l, xi_u, RESign, this};
        gsl_integration_qag(&int_gsl_func,p_l,p_u,epsabs,epsrel,lim,QAG_KEY,gsl_ws_CH,&avaCH_FVA_passing_test, &error);

        real_t p_i  =  (p_l + p_u)  / 2.;
        real_t xi_j = (xi_l + xi_u) / 2.;

        real_t Delta_p  = p_u - p_l;
        real_t Delta_xi = xi_u - xi_l;

        avaCH_FVA_passing_test *= VpVol / (Vp * Delta_p * Delta_xi);

        printf("\nBA=%.4e, FVA=%.4e, diff=%.4e\n", avaCH_BA, avaCH_FVA_passing_test, abs(avaCH_BA_passing_test - avaCH_FVA_passing_test));
    } 
    */
    
    return avaCH_BA;
}