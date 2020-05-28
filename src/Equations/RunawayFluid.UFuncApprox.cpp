

// For GSL functions: partial contributions to evaluateUAtP
struct UFuncParams {len_t ir; real_t A; FVM::RadialGrid *rGrid;};
real_t UFrictionTermIntegrand(real_t xi0, void *par){
    struct UFuncParams *params = (struct UFuncParams *) par;
    len_t ir = params->ir;
    real_t A = params->A;
    FVM::RadialGrid *rGrid = params->rGrid;
    if(xi0==0)
        return exp(-A);
    std::function<real_t(real_t,real_t,real_t)> FrictionTermFunc = [xi0](real_t BOverBmin, real_t , real_t )
                            {return sqrt(xi0*xi0/(1-BOverBmin*(1-xi0*xi0)))*BOverBmin;};
    return rGrid->CalculateFluxSurfaceAverage(ir,false, FrictionTermFunc)*exp(-A*(1-xi0));
}

real_t USynchrotronTermIntegrand(real_t xi0, void *par){
    struct UFuncParams *params = (struct UFuncParams *) par;
    len_t ir = params->ir;
    real_t A = params->A;
    FVM::RadialGrid *rGrid = params->rGrid;
    if(xi0==0)
        return exp(-A);
    std::function<real_t(real_t,real_t,real_t)> SynchrotronTermFunc = [xi0](real_t BOverBmin, real_t , real_t )
                            {return (1-xi0*xi0)*sqrt(xi0*xi0/(1-BOverBmin*(1-xi0*xi0))) *BOverBmin*BOverBmin*BOverBmin*BOverBmin;};
    return rGrid->CalculateFluxSurfaceAverage(ir,false, SynchrotronTermFunc)*exp(-A*(1-xi0));
}


// Evaluates the effective momentum flow U accounting for electric field, collisional friction and radiation reaction 
real_t RunawayFluid::evaluateApproximateUAtP(real_t p, void *par){
    struct UContributionParams *params = (struct UContributionParams *) par;
    struct CollisionQuantity::collqty_settings *collSettingsForEc = params->collSettingsForEc;
    bool rFluxGrid = params->rFluxGrid;
    real_t Eterm = params->Eterm;
    len_t ir = params->ir;
    gsl_integration_workspace *gsl_ad_w = params->gsl_ad_w;
    FVM::RadialGrid *rGrid = params->rGrid;
    PitchScatterFrequency *nuD = params->nuD;
    SlowingDownFrequency *nuS = params->nuS;

    real_t Bmin,Bmax;
    if(rFluxGrid){
        Bmin = rGrid->GetBmin_f(ir);
        Bmax = rGrid->GetBmax_f(ir);    
    }else{
        Bmin = rGrid->GetBmin(ir);
        Bmax = rGrid->GetBmax(ir);
    }

    const real_t B2avgOverBmin2 = rGrid->GetFSA_B2(ir);
    
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(B2avgOverBmin2); 
    real_t xiT = sqrt(1-Bmin/Bmax);
//    real_t xiT = -1;
    real_t pNuD = p*nuD->evaluateAtP(ir,p,collSettingsForEc);
    real_t A = 2*E/pNuD;
    real_t Econtrib;
    if(A==0)
        Econtrib = 0.5*E*(1-xiT*xiT);
    else 
        Econtrib = E/(A*A) *( A-1 - exp(-A*(1-xiT))*(A*xiT -1) );

    real_t FrictionTerm = p*nuS->evaluateAtP(ir,p,collSettingsForEc);
    
    UFuncParams FuncParams = {ir, A, rGrid};
    gsl_function UIntegrandFunc;

    UIntegrandFunc.function = &(UFrictionTermIntegrand);
    UIntegrandFunc.params = &FuncParams;
    real_t abserr;
    real_t frictionIntegral;
//    real_t frictionIntegralPredict = (1-exp(-A))/A;

    real_t epsabs=0, epsrel=1e-3, lim = gsl_ad_w->limit;
    gsl_integration_qags(&UIntegrandFunc, xiT,1.0,epsabs,epsrel,lim,gsl_ad_w, &frictionIntegral, &abserr);
    real_t FrictionContrib = -FrictionTerm * frictionIntegral;

    real_t SynchrotronTerm = p*sqrt(1+p*p)* Constants::ec * Constants::ec * Constants::ec * Constants::ec * Bmin * Bmin
                            / ( 6 * M_PI * Constants::eps0 * Constants::me * Constants::me * Constants::me
                                * Constants::c * Constants::c * Constants::c); 
                               
    UIntegrandFunc.function = &(USynchrotronTermIntegrand);
    real_t synchrotronIntegral;
//    real_t synchrotronIntegralPredict = exp(-A)*(-A*A+2*(A-1)*exp(A)+2)/(A*A*A);
    gsl_integration_qags(&UIntegrandFunc, xiT,1.0,epsabs,epsrel,lim,gsl_ad_w, &synchrotronIntegral, &abserr);
    real_t SynchrotronContrib = -SynchrotronTerm * synchrotronIntegral;
    
   return -(Econtrib + FrictionContrib + SynchrotronContrib) / frictionIntegral;
}



real_t RunawayFluid::evaluateApproximatePitchDistribution(len_t ir, real_t xi0, real_t p, real_t Eterm, CollisionQuantity::collqty_settings *inSettings){

    const real_t B2avgOverBmin2 = rGrid->GetFSA_B2(ir);
    real_t E = Constants::ec * Eterm / (Constants::me * Constants::c) * sqrt(B2avgOverBmin2); 
    real_t pNuD = p*nuD->evaluateAtP(ir,p,inSettings);
    real_t A = 2*E/pNuD;

    return exp(-A*(1-xi0));
}
        