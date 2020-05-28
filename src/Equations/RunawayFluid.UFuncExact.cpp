/**
 * U(r,p) is the pitch-angle averaged momentum-advection coefficient,
 * U = < {A^p}f > / <f>,
 * where <...> = int ... d(xi0).
 * The calculation is documented under doc/notes/theory.pdf in 
 * Section 2 (under the heading 'Bounce-averaged effective field') 
 */
struct distExponentParams {len_t ir; FVM::RadialGrid *rGrid;};
real_t distExponentIntegral(real_t xi0, void *par){
    struct distExponentParams *params = (struct distExponentParams *) par;
    std::function<real_t(real_t,real_t,real_t)> xiFunc = [xi0](real_t BOverBmin, real_t , real_t )
                            {return sqrt(1 - BOverBmin*(1-xi0*xi0));};
    len_t ir = params->ir;
    FVM::RadialGrid *rGrid = params->rGrid;
    real_t signXi0 = ( (xi0>0) - (xi0<0));
    real_t xiAvg = signXi0 * rGrid->CalculateFluxSurfaceAverage(ir,false, xiFunc);
    return xi0/xiAvg;
}



real_t RunawayFluid::evaluateAnalyticPitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, gsl_integration_workspace *gsl_ad_w){
    const real_t Bmin = rGrid->GetBmin(ir);
    const real_t Bmax = rGrid->GetBmax(ir);
    const real_t B2avgOverBmin2 = rGrid->GetFSA_B2(ir);
    real_t xiT = sqrt(1-Bmin/Bmax);
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(B2avgOverBmin2); 

//    const CollisionQuantity::collqty_settings *collQtySettings = rf->GetSettings();
    real_t pNuD = p*nuD->evaluateAtP(ir,p,collQtySettings->collfreq_type,OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL);    
    real_t A = 2*E/pNuD;

    // This block carries defines the integration int(xi0/<xi(xi0)> dxi0, xi1, x2) 
    //////////////////////////////
    gsl_function GSL_func;
    distExponentParams params = {ir,rGrid};
    GSL_func.function = &(distExponentIntegral);
    GSL_func.params = &params;
    real_t abserr;
    real_t epsabs = 0, epsrel = 1e-4, lim = gsl_ad_w->limit;
    #define F(xi1,xi2,pitchDist) gsl_integration_qags(&GSL_func, xi1,xi2,epsabs,epsrel,lim,gsl_ad_w, &pitchDist, &abserr)
    //////////////////////////////    

    real_t dist1 = 0;
    real_t dist2 = 0;

    real_t thresholdToNeglectTrappedContribution = 1e-6;
    if ( (xi0>xiT) || (xiT<thresholdToNeglectTrappedContribution) )
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



/*
The function takes a xi0 and a lambda expression Func (and other needed helper parameters) and 
returns the contribution to the integrand in the U function, i.e. V'{Func}*exp(-...),
where exp(-...)(xi0) is the analytical pitch-angle distribution, and V'{Func} the 
bounce integral of Func.
*/

real_t UPartialContribution(real_t xi0, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    
    RunawayFluid *rf = params->rf; 
    FVM::RadialGrid *rGrid = params->rGrid; 
    len_t ir = params->ir;
    real_t p = params->p;
    bool rFluxGrid = params->rFluxGrid;
    gsl_integration_workspace *gsl_ad_w = params->gsl_ad_w;
    real_t E = params->Eterm;
    std::function<real_t(real_t,real_t)> BAFunc = [xi0,params](real_t x,real_t y){return params->Func(xi0,x,y);};
    
    return rGrid->evaluatePXiBounceIntegralAtP(ir,p,xi0,rFluxGrid,BAFunc,gsl_ad_w)
        * rf->evaluateAnalyticPitchDistribution(ir,xi0,p,E,gsl_ad_w);    
}

real_t RunawayFluid::evaluateNegUAtP(real_t p, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    params->p = p;
    FVM::RadialGrid *rGrid = params->rGrid;
    len_t ir = params->ir;
    bool rFluxGrid = params->rFluxGrid;
    real_t Eterm = params->Eterm;
    gsl_integration_workspace *gsl_ad_w = params->gsl_ad_w;
    SlowingDownFrequency *nuS = params->nuS;
    RunawayFluid *rf = params->rf;

    real_t Bmin,Bmax;
    if(rFluxGrid){
        Bmin = rGrid->GetBmin_f(ir);
        Bmax = rGrid->GetBmax_f(ir);    
    }else{
        Bmin = rGrid->GetBmin(ir);
        Bmax = rGrid->GetBmax(ir);
    }
    const real_t sqrtB2avgOverBavg = sqrt(rGrid->GetFSA_B2(ir)) / rGrid->GetFSA_B(ir);
    real_t xiT = sqrt(1-Bmin/Bmax);
    if(xiT < 1e-6)
        xiT = 0;

    // Evaluates the contribution from electric field term A^p coefficient
    std::function<real_t(real_t,real_t,real_t)> FuncElectric = 
            [](real_t xi0, real_t /*BOverBmin*/, real_t xiOverXi0 ){return xi0*xiOverXi0;};

    params->Func = FuncElectric;
    real_t EContrib,EContrib1,EContrib2,error;
    real_t Efactor = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrtB2avgOverBavg; 
    real_t epsabs = 0, epsrel = 5e-3, lim = gsl_ad_w->limit; 
    gsl_function GSL_func;
    GSL_func.function = &(UPartialContribution);
    GSL_func.params = params;
    if(xiT){
        gsl_integration_qags(&GSL_func,-1,-xiT,epsabs,epsrel,lim,gsl_ad_w,&EContrib1,&error);
        gsl_integration_qags(&GSL_func,xiT,1,epsabs,epsrel,lim,gsl_ad_w,&EContrib2,&error);
        EContrib = EContrib1 + EContrib2;
    }else
        gsl_integration_qags(&GSL_func,-1,1,epsabs,epsrel,lim,gsl_ad_w,&EContrib,&error);
    EContrib *= Efactor;

    // Evaluates the contribution from slowing down term A^p coefficient
    std::function<real_t(real_t,real_t,real_t)> FuncUnity = 
            [](real_t,real_t,real_t){return 1;};
    params->Func = FuncUnity;    
    real_t UnityContrib,UnityContrib1,UnityContrib2,UnityContrib3;
    if(xiT){
        gsl_integration_qags(&GSL_func,-1,-xiT,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib1,&error);
        gsl_integration_qags(&GSL_func,-xiT,xiT,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib2,&error);
        gsl_integration_qags(&GSL_func,xiT,1,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib3,&error);
        UnityContrib = UnityContrib1 + UnityContrib2 + UnityContrib3;
    } else 
        gsl_integration_qags(&GSL_func,-1,1,epsabs,epsrel,lim,gsl_ad_w,&UnityContrib,&error);

    const CollisionQuantity::collqty_settings *collQtySettings = rf->GetSettings();
    real_t NuSContrib = -p*nuS->evaluateAtP(ir,p,collQtySettings->collfreq_type,OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL)* UnityContrib;


    // Evaluates the contribution from synchrotron term A^p coefficient
    std::function<real_t(real_t,real_t,real_t)> FuncSynchrotron = 
            [](real_t xi0, real_t BOverBmin, real_t){return (1-xi0*xi0)*BOverBmin*BOverBmin*BOverBmin;};
    params->Func = FuncSynchrotron;
    real_t SynchrotronFactor = -p*sqrt(1+p*p)* Constants::ec * Constants::ec * Constants::ec * Constants::ec * Bmin * Bmin
                            / ( 6 * M_PI * Constants::eps0 * Constants::me * Constants::me * Constants::me
                                * Constants::c * Constants::c * Constants::c); 

    real_t SynchContrib,SynchContrib1,SynchContrib2,SynchContrib3;
    if(xiT){
        gsl_integration_qags(&GSL_func,-1,-xiT,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib1,&error);
        gsl_integration_qags(&GSL_func,-xiT,xiT,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib2,&error);
        gsl_integration_qags(&GSL_func,xiT,1,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib3,&error);
        SynchContrib = SynchContrib1 + SynchContrib2 + SynchContrib3;
    } else 
        gsl_integration_qags(&GSL_func,-1,1,epsabs,epsrel,lim,gsl_ad_w,&SynchContrib,&error);

    SynchContrib *= SynchrotronFactor; 

    return -(EContrib + NuSContrib + SynchContrib) / UnityContrib;

}
