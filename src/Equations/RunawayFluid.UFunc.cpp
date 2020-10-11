/**
 * U(r,p) is the pitch-angle averaged momentum-advection coefficient,
 * U = < {A^p}f > / <f>,
 * where <...> = int ... d(xi0).
 * The calculation is documented under doc/notes/theory.pdf in 
 * Section 2 (under the heading 'Bounce-averaged effective field') 
 */

const real_t THRESHOLD_NEGLECT_TRAPPING = 1e-6;

/**
 * Returns xi0/<xi> (the integral of which appears in AnalyticPitchDistribution).
 */
struct distExponentParams {len_t ir; FVM::RadialGrid *rGrid;};
real_t distExponentIntegral(real_t xi0, void *par){
    struct distExponentParams *params = (struct distExponentParams *) par;
    std::function<real_t(real_t,real_t,real_t)> xiFunc = [xi0](real_t BOverBmin, real_t , real_t )
                            {return sqrt(1 - BOverBmin*(1-xi0*xi0));};
    len_t ir = params->ir;
    FVM::RadialGrid *rGrid = params->rGrid;
    real_t signXi0 = ( (xi0>0) - (xi0<0));
    real_t xiAvg = signXi0 * rGrid->CalculateFluxSurfaceAverage(ir,FVM::FLUXGRIDTYPE_DISTRIBUTION, xiFunc);
    return xi0/xiAvg;
}


/**
 * Calculates the (semi-)analytic pitch-angle distribution predicted in the 
 * near-threshold regime, where the momentum flux is small compared 
 * to the characteristic pitch flux, and we obtain the approximate 
 * kinetic equation phi_xi = 0.
 */
real_t RunawayFluid::evaluateAnalyticPitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings, gsl_integration_workspace *gsl_ad_w){
    const real_t Bmin = rGrid->GetBmin(ir);
    const real_t Bmax = rGrid->GetBmax(ir);
    const real_t B2avgOverBmin2 = rGrid->GetFSA_B2(ir);
    real_t xiT = sqrt(1-Bmin/Bmax);
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(B2avgOverBmin2); 

    real_t pNuD = p*nuD->evaluateAtP(ir,p,inSettings);    
    real_t A = 2*E/pNuD;

    // This block carries defines the integration int(xi0/<xi(xi0)> dxi0, xi1, x2) 
    //////////////////////////////
    gsl_function GSL_func;
    distExponentParams params = {ir,rGrid};
    GSL_func.function = &(distExponentIntegral);
    GSL_func.params = &params;
    real_t abserr;
    real_t epsabs = 0, epsrel = 3e-3, lim = gsl_ad_w->limit;
    #define F(xi1,xi2,pitchDist) gsl_integration_qag(&GSL_func, xi1,xi2,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w, &pitchDist, &abserr)
    //////////////////////////////    

    real_t dist1 = 0;
    real_t dist2 = 0;

    if ( (xi0>xiT) || (xiT<THRESHOLD_NEGLECT_TRAPPING) )
        F(xi0,1.0,dist1);
    else if ( (-xiT <= xi0) && (xi0 <= xiT) )
        F(xiT,1.0,dist1);
    else{ // (xi0 < -xiT)
        F(xi0,-xiT,dist1);
        F(xiT,1.0,dist2);
    }
    
    #undef F
    
    return exp(-A*(dist1+dist2));
}

real_t RunawayFluid::evaluatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings, gsl_integration_workspace *gsl_ad_w){
    if(Eceff_mode == OptionConstants::COLLQTY_ECEFF_MODE_SIMPLE)
        return evaluateApproximatePitchDistribution(ir,xi0,p,Eterm,inSettings);
    else
        return evaluateAnalyticPitchDistribution(ir,xi0,p,Eterm,inSettings,gsl_ad_w);
}

/*
The function takes a xi0 and a lambda expression Func (and other needed helper parameters) and 
returns the contribution to the integrand in the U function, i.e. V'{Func}*exp(-...),
where exp(-...)(xi0) is the analytical pitch-angle distribution, and V'{Func} the 
bounce integral of Func.
*/

real_t UPartialContribution(real_t xi0, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    CollisionQuantity::collqty_settings *collSettingsForEc = params->collSettingsForEc;
    RunawayFluid *rf = params->rf; 
    FVM::RadialGrid *rGrid = params->rGrid; 
    len_t ir = params->ir;
    real_t p = params->p;
    FVM::fluxGridType fluxGridType = params->fgType;
    gsl_integration_workspace *gsl_ad_w = params->gsl_ad_w;
    real_t E = params->Eterm;
    std::function<real_t(real_t,real_t,real_t,real_t)> BAFunc = 
        [xi0,params](real_t xiOverXi0,real_t BOverBmin,real_t /*ROverR0*/,real_t /*NablaR2*/)
            {return params->Func(xi0,BOverBmin,xiOverXi0);};
    
    return rGrid->EvaluatePXiBounceIntegralAtP(ir,p,xi0,fluxGridType,BAFunc)
        * rf->evaluatePitchDistribution(ir,xi0,p,E,collSettingsForEc, gsl_ad_w);    
}

/**
 * Evaluates -U(p) at given Eterm.
 */
real_t RunawayFluid::UAtPFunc(real_t p, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    params->p = p;
    FVM::RadialGrid *rGrid = params->rGrid;
    len_t ir = params->ir;
    FVM::fluxGridType fluxGridType = params->fgType;
    real_t Eterm = params->Eterm;
    gsl_integration_workspace *gsl_ad_w = params->gsl_ad_w;
    SlowingDownFrequency *nuS = params->nuS;
    CollisionQuantity::collqty_settings *collSettingsForEc = params->collSettingsForEc;
    int QAG_KEY = params->QAG_KEY;
    real_t Bmin,Bmax;
    if(fluxGridType == FVM::FLUXGRIDTYPE_RADIAL){
        Bmin = rGrid->GetBmin_f(ir);
        Bmax = rGrid->GetBmax_f(ir);    
    }else{
        Bmin = rGrid->GetBmin(ir);
        Bmax = rGrid->GetBmax(ir);
    }
    const real_t sqrtB2avgOverBavg = sqrt(rGrid->GetFSA_B2(ir)) / rGrid->GetFSA_B(ir);
    real_t xiT = sqrt(1-Bmin/Bmax);
    if(xiT < THRESHOLD_NEGLECT_TRAPPING)
        xiT = 0;

    // Evaluates the contribution from electric field term A^p coefficient
    std::function<real_t(real_t,real_t,real_t)> FuncElectric = 
            [](real_t xi0, real_t /*BOverBmin*/, real_t xiOverXi0 ){return xi0*xiOverXi0;};

    params->Func = FuncElectric;
    real_t EContrib, error;
    real_t Efactor = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrtB2avgOverBavg; 
    real_t epsabs = 0, epsrel = 1e-3, lim = gsl_ad_w->limit; 
    gsl_function GSL_func;
    GSL_func.function = &(UPartialContribution);
    GSL_func.params = params;
    if(xiT){
        real_t EContrib1, EContrib2;
        gsl_integration_qag(&GSL_func,-1,-xiT,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&EContrib1,&error);
        gsl_integration_qag(&GSL_func,xiT,1,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&EContrib2,&error);
        EContrib = EContrib1 + EContrib2;
    }else
        gsl_integration_qag(&GSL_func,-1,1,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&EContrib,&error);
    EContrib *= Efactor;

    // Evaluates the contribution from slowing down term A^p coefficient
    std::function<real_t(real_t,real_t,real_t)> FuncUnity = 
            [](real_t,real_t,real_t){return 1;};
    params->Func = FuncUnity;    
    real_t UnityContrib;
    if(xiT){
        real_t UnityContrib1, UnityContrib2, UnityContrib3;
        gsl_integration_qag(&GSL_func,-1,-xiT,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&UnityContrib1,&error);
        gsl_integration_qag(&GSL_func,0,xiT,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&UnityContrib2,&error);
        gsl_integration_qag(&GSL_func,xiT,1,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&UnityContrib3,&error);
        UnityContrib = UnityContrib1 + UnityContrib2 + UnityContrib3;
    } else 
        gsl_integration_qag(&GSL_func,-1,1,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&UnityContrib,&error);

    real_t NuSContrib = -p*nuS->evaluateAtP(ir,p,collSettingsForEc) * UnityContrib;


    // Evaluates the contribution from synchrotron term A^p coefficient
    std::function<real_t(real_t,real_t,real_t)> FuncSynchrotron = 
            [](real_t xi0, real_t BOverBmin, real_t){return (1-xi0*xi0)*BOverBmin*BOverBmin*BOverBmin;};
    params->Func = FuncSynchrotron;
    real_t SynchrotronFactor = -p*sqrt(1+p*p)* Constants::ec * Constants::ec * Constants::ec * Constants::ec * Bmin * Bmin
                            / ( 6 * M_PI * Constants::eps0 * Constants::me * Constants::me * Constants::me
                                * Constants::c * Constants::c * Constants::c); 

    real_t SynchContrib;
    if(xiT){
        real_t SynchContrib1, SynchContrib2, SynchContrib3;
        gsl_integration_qag(&GSL_func,-1,-xiT,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&SynchContrib1,&error);
        gsl_integration_qag(&GSL_func,0,xiT,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&SynchContrib2,&error);
        gsl_integration_qag(&GSL_func,xiT,1,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&SynchContrib3,&error);
        SynchContrib = SynchContrib1 + SynchContrib2 + SynchContrib3;
    } else 
        gsl_integration_qag(&GSL_func,-1,1,epsabs,epsrel,lim,QAG_KEY,gsl_ad_w,&SynchContrib,&error);

    SynchContrib *= SynchrotronFactor; 

    return -(EContrib + NuSContrib + SynchContrib) / UnityContrib;

}













/**
 * Same as evaluteAnalyticPitchDistribution, but approximating
 * xi0/<xi> = 1 for passing and 0 for trapped (thus avoiding the 
 * need for the numerical integration).
 */
real_t RunawayFluid::evaluateApproximatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings){
    const real_t Bmin = rGrid->GetBmin(ir);
    const real_t Bmax = rGrid->GetBmax(ir);
    const real_t B2avgOverBmin2 = rGrid->GetFSA_B2(ir);
    real_t xiT = sqrt(1-Bmin/Bmax);
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(B2avgOverBmin2); 

//    const CollisionQuantity::collqty_settings *collQtySettings = rf->GetSettings();
    real_t pNuD = p*nuD->evaluateAtP(ir,p,inSettings);    
    real_t A = 2*E/pNuD;

    real_t dist1 = 0;
    real_t dist2 = 0;

    if ( (xi0>xiT) || (xiT<THRESHOLD_NEGLECT_TRAPPING) )
        dist1 = 1-xi0;
    else if ( (-xiT <= xi0) && (xi0 <= xiT) )
        dist1 = 1-xiT;
    else{ // (xi0 < -xiT)
        dist1 = 1-xiT;
        dist2 = -xiT - xi0;
    }
        
    return exp(-A*(dist1+dist2));
}