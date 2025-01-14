/**
 * In this file, methods used for the calculation of the bounce averaged
 * Chiu-Harvey avalanche source are collected.
 */

/**
 * Helper function for avalanche bounce average calculation: 
 * returns the function 
 * xi0_min = sqrt((gamma - 1) / (gamma + 1) 
 *                          * (gamma_{in,max} + 1) / (gamma_{in,max} - 1))
 */
real_t ximin(real_t gamma, real_t gamma_max){
    return sqrt((gamma - 1) / (gamma + 1) * (gamma_max + 1) / (gamma_max - 1));
}

/**
 * Helper function for avalanche bounce average calculation: 
 * returns the function 
 * xi0_max = sqrt(gamma / (gamma + 1))
 */
real_t ximax(real_t gamma){
    return sqrt(gamma / (gamma + 1));
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
    
    if (ximax(gamma) < ximin(gamma, gamma_max))
        return 0;
    if (ximax(gamma) < xiOverXi0*abs(xi0))
        return 0;
    if (ximin(gamma, gamma_max) > xiOverXi0*abs(xi0))
        return 0;
    
    real_t xi2 = 1. - (1. - xi0*xi0) * BOverBmin;
    real_t denom = xi2 - (gamma - 1.) / (gamma + 1.);
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
    int_t iTerm = params->iTerm;
    FluxSurfaceAverager *FSA = params->FSA;
    
    if(xi0 < ximin(gamma, gamma_max))
        return 0;
    
    avParams avg_params[5] = {xi0, gamma, gamma_max, iTerm, BminOverBmax};
    
    bool doQAGS = true;
    
    real_t factor_BA = FSA->CalculatePXiBounceAverageAtP(ir, xi0, fgt, &BA_CH, avg_params, nullptr, doQAGS);
    
    if ((1-xi0*xi0) > BminOverBmax)
        factor_BA *= 0.5;
    
    if (iTerm < 3){
        return xi0*xi0 / ((1 - xi0*xi0) * (1 - xi0*xi0)) * factor_BA;
    }

    return xi0*xi0 / (1 - xi0*xi0) * factor_BA;
}

/*
 * Integrand from discretization in momentum.
 */
real_t FluxSurfaceAverager::integrandP(real_t p, real_t gamma_max, len_t ir, real_t xi_l, real_t xi_u, real_t BminOverBmax, fluxGridType fgt){
    real_t gamma = sqrt(1+p*p);
    real_t epsabs = 0, epsrel = 1e-3, error;
    
    gsl_function int_gsl_func;
    int_gsl_func.function = &(integrandXi);

    real_t terms[6]; 
    for(int_t i=0; i<6; i++){
        intXiParams intXi_params = {ir, fgt, gamma, gamma_max, BminOverBmax, i, this};
        int_gsl_func.params = &intXi_params;
        gsl_integration_qags(&int_gsl_func,xi_l,xi_u,epsabs,epsrel, gsl_ws_CH->limit,gsl_ws_CH,&terms[i], &error);
        //len_t nevals;
        //gsl_integration_cquad(&int_gsl_func,xi_l,xi_u,epsabs,epsrel,gsl_ws_CH,&terms[i],&error,&nevals);
        
    }



    return 1/(gamma*p)*
            (4*(gamma-1)/((gamma+1)*(gamma+1)*(gamma+1)*(gamma+1)*(gamma+1))*terms[0]
              +8/((gamma+1)*(gamma+1)*(gamma+1)*(gamma+1))*terms[1]
              +4/((gamma-1)*(gamma+1)*(gamma+1)*(gamma+1))*terms[2]
              -1*(gamma-1)*(gamma*gamma-2)/((gamma+1)*(gamma+1)*(gamma+1))*terms[3]
              +2*(gamma*gamma-gamma-3)/((gamma+1)*(gamma+1))*terms[4]
              -(gamma*gamma-2*gamma+4)/((gamma-1)*(gamma+1))*terms[5]);
}
/**
 * Returns the bounce and cell averaged Chiu-Harvey avalanche
 * source term.
 * 
 * Parameters:
 *       ir: radial grid point
 *      p_l: momentum on the lower cell face
 *      p_u: momentum on the upper cell face
 *    p_max: maxmimum momentum of RE grid
 *     xi_l: xi on the lower cell face
 *     xi_u: xi on the upper cell face
 *   RESign: sign of xi of the incident REs (+1 or -1).
 *           Is used to flip the pitch of the source
 *      fgt: fluxGridType object needed for bounce averaging.
 */
real_t FluxSurfaceAverager::EvaluateAvalancheCHBounceAverage(len_t ir, real_t p_i, real_t p_max, real_t xi_l, real_t xi_u,  fluxGridType fgt){
    real_t Bmin = GetBmin(ir, FLUXGRIDTYPE_DISTRIBUTION,nullptr);
    real_t Bmax = GetBmax(ir, FLUXGRIDTYPE_DISTRIBUTION,nullptr);
    real_t BminOverBmax;
    real_t xi_l_min = xi_l; 
    real_t Delta_xi = xi_u - xi_l;
    real_t xi_j = (xi_u + xi_l) / 2.;
    real_t gamma_i = sqrt(1+p_i*p_i);
    real_t gamma_max = sqrt(1+p_max*p_max);
    real_t xi_min = ximin(gamma_i, gamma_max);
    real_t xi_max = ximax(gamma_i);
    
    if (xi_j < 0){
        xi_j *= -1;
        real_t xi_l_temp = xi_l;
        xi_l = -xi_u;
        xi_u = -xi_l_temp;
    }
    
    if(Bmin==Bmax)
        BminOverBmax=1;
    else {
        BminOverBmax = Bmin/Bmax;
        real_t xi_trapped = sqrt(1 - BminOverBmax);
        if (xi_trapped < xi_l)
            xi_l_min = sqrt(1 - (1 - xi_l*xi_l)/BminOverBmax);
        else 
            xi_l_min = 0;
    }

    if( xi_max <= xi_l_min )
        return 0;
    else if( xi_min >= xi_u )
        return 0;
    else if ( xi_max < xi_min)
        return 0;
    
    return 1 / (abs(xi_j) * Delta_xi) * integrandP(p_i, gamma_max, ir, xi_l, xi_u, BminOverBmax, fgt);
}