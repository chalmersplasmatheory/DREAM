/**
 * Implementation of equation term representing the runaway generation rate
 * due to hottail when using an analytic distribution function. The generation
 * rate is calculated according to section 4.1 of Ida Svenningsson's MSc thesis
 * (https://hdl.handle.net/20.500.12380/300899).
 */

#include "DREAM/Equations/Fluid/HottailRateTermLowZ.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
HottailRateTermLowZ::HottailRateTermLowZ(           
    FVM::Grid *grid, AnalyticDistributionHottail *dist, FVM::UnknownQuantityHandler *unknowns,
    IonHandler *ionHandler, CoulombLogarithm *lnL, RunawayFluid *runawayFluid, const real_t* T_final, real_t sf
) : HottailRateTerm(grid, dist, unknowns,sf), lnL(lnL), ionHandler(ionHandler), runawayFluid(runawayFluid), T_final(T_final), 
    id_ncold(unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD)),
    id_Efield(unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD)), // Used for constant, not in jacobian.
    id_tau(unknowns->GetUnknownID(OptionConstants::UQTY_TAU_COLL)),
    id_Tcold(unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD)), // Needed for derivatives in jacobian
    id_johm(unknowns->GetUnknownID(OptionConstants::UQTY_J_OHM)), // Used for constant, not in jacobian. 
    id_ni(unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES)) // Needed for derivatives in jacobian
{
    SetName("HottailRateTermLowZ");

    AddUnknownForJacobian(unknowns,id_ncold);
    AddUnknownForJacobian(unknowns,id_tau);
    AddUnknownForJacobian(unknowns,id_ni); //Also affects Zeff, now we neglect this contribution. 
    AddUnknownForJacobian(unknowns,id_Tcold);
    
    this->w = gsl_integration_workspace_alloc(nGslIntervals);
    this->rGrid = grid->GetRadialGrid();
    
    this->GridRebuilt();
}

/**
 * Destructor
 */
HottailRateTermLowZ::~HottailRateTermLowZ(){
    Deallocate();
    gsl_integration_workspace_free(w);
}           

/**
 * Called after the grid is rebuilt; (re)allocates memory
 * for all quantities
 */
bool HottailRateTermLowZ::GridRebuilt(){
    this->HottailRateTerm::GridRebuilt();

    Deallocate();
            
    this->params = new ParamStruct();
    params->rGrid = this->rGrid; 
    params->dist = distHT;
    
    Epar_prev = new real_t[nr];
    Epar = new real_t[nr]; 
    
    params->sp_cond = new real_t[nr];
    params->Zeff = new real_t[nr];
    params->ncold = new real_t[nr];
    params->tau = new real_t[nr];
    params->j0 = new real_t[nr];
    params->lnL = new real_t[nr];
    params->Ec = new real_t[nr];
    
    integraloff0inGamma = new real_t[nr];

    return true;
}

/**
 * Rebuilds quantities used by this equation term.
 * Note that the equation term uses the _energy distribution_
 * from AnalyticDistributionHottail which differs from the 
 * full distribution function by a factor of 4*pi
 */
void HottailRateTermLowZ::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler* u) {
    if (Epar == nullptr){
    	Epar = new real_t[nr];
    	const real_t *Einit = u->GetUnknownInitialData(id_Efield);
    	for(len_t ir=0; ir<nr; ir++){
    	    Epar[ir] = Einit[ir];
    	}
    }
    this->dt = dt;
    bool newTimeStep = (t!=tPrev);
    if(newTimeStep)
        tPrev = t;
    
    params->Zeff = ionHandler->GetZeff(); 
    params->ncold = u->GetUnknownData(id_ncold);
    params->tau = u->GetUnknownData(id_tau);
    params->j0 = u->GetUnknownInitialData(id_johm);
    
    for(len_t ir=0; ir<nr; ir++){
    
        params->ir = ir;
        params->lnL[ir] = lnL->evaluateAtP(ir,0); 
        params->Ec[ir] = 4*M_PI*(params->ncold[ir])*(params->lnL[ir])*Constants::r0*Constants::r0*Constants::c*Constants::c*Constants::me / Constants::ec; 
        params->sp_cond[ir] = runawayFluid->evaluateBraamsElectricConductivity(ir, T_final[ir], params->Zeff[ir]);
        
        if(newTimeStep)
            Epar_prev[ir] = Epar[ir];
        Epar[ir] = evaluate_Epar(params);
        real_t dotEpar = (Epar[ir] - Epar_prev[ir]) / dt;
         if (dotEpar < 0){ // ensure non-negative runaway rate
            dotEpar = 0;
            }
        real_t integral_f0 = integralf0(Epar[ir], params);
        integraloff0inGamma[ir] = integral_f0; // Save for jacobian to not need to calc again
        gamma[ir] = dotEpar* integral_f0* (params->Ec[ir]) / (Epar[ir]*Epar[ir]); // factor 4pi built into def. of f0 --> no 4pi here
    }
}


/**
* Evaluates f_0(t,p) (4pi times the f0 def. in Ida Svenningssons MSc thesis) 
*/
real_t HottailRateTermLowZ::evaluate_f0(real_t p, void *par){ 
    
    ParamStruct *pars = (ParamStruct*)par;
    len_t ir = pars->ir;
    real_t tau = pars->tau[ir];
    
    real_t f0 = pars->dist->evaluateEnergyDistributionFromTau(ir,p,tau); 
    return f0; 
}


/**
* Evaluates other part of integrand for E_parallell at momentum p
*/
real_t HottailRateTermLowZ::partialIntegrandForEpar(real_t p){
    real_t p2 = p*p; 
    // p^5 = p2*p2*p
    // p^7 = p2*p2*p2*p
    real_t partialIntegrand = (6*p2*p2*p + 4*p2*p2*p2*p) / ((1+p2)*(1+p2));
    return partialIntegrand;
}

/**
* Evaluates total integrand for E_parallell at momentum p
*/
real_t HottailRateTermLowZ::totalIntegrandForEpar(real_t p, void *par){
    real_t totalIntegrand = evaluate_f0(p, par)*partialIntegrandForEpar(p);
    return totalIntegrand;
}

/**
* Evaluates integral in expression E_parallell at momentum p
*/
real_t HottailRateTermLowZ::integralEpar(struct ParamStruct * intparams){ 

    real_t result, error;
    
    gsl_function Func;
    Func.function = &totalIntegrandForEpar;
    Func.params = intparams;
    
    gsl_integration_qagiu(&Func, 0, ABSTOL_FOR_INT, RELTOL_FOR_INT, nGslIntervals, w, &result, &error);
    
    return result;
}

/**
* Calculates E_parallell
*/
real_t HottailRateTermLowZ::evaluate_Epar(struct ParamStruct * intparams){

    len_t ir = intparams->ir;
    real_t Ec = intparams->Ec[ir];
    real_t sp_cond = intparams->sp_cond[ir];
    real_t j0 = intparams->j0[ir];
    real_t Zeff = intparams->Zeff[ir];
    
    real_t integ = integralEpar(intparams);
    real_t Epar_ir = Ec*j0 / (sp_cond*Ec + integ*Constants::ec*Constants::c / (3*(1+Zeff))); // Removed 4pi factor in second term in denominator since included in f0
    
    return Epar_ir;
}

/**
* Evaluates integral of f0 in expression for gamma, includes the factor 4pi
*/
real_t HottailRateTermLowZ::integralf0(real_t Epar_ir, struct ParamStruct * intparams){ 
    
    real_t result, error;
    len_t ir = intparams->ir;
    real_t lowerlimit = sqrt((intparams->Ec[ir]) / Epar_ir);
    
    gsl_function Func;
    Func.function = &evaluate_f0;
    Func.params = intparams;
    
    gsl_integration_qagiu(&Func, lowerlimit, ABSTOL_FOR_INT, RELTOL_FOR_INT, nGslIntervals, w, &result, &error);
    
    return result;
}

/**
* Return E-field
*/
const real_t* HottailRateTermLowZ::GetElectricField(){
    return this->Epar;
}


/**************************************************************/
/****************** Functions for Jacobian ********************/
/**************************************************************/
// Also use some of the fcns above


/**
* Derivative of Ec wrt unknowns
*/
real_t HottailRateTermLowZ::dEcdu(struct dGammadu_Params * pars){

    real_t prefactor = 4*M_PI*Constants::r0*Constants::r0*Constants::c*Constants::c*Constants::me / Constants::ec;
    real_t var;
    
    if (pars->id_unknown == id_ncold) {
        var = pars->lnL;
    } else if (pars->id_unknown == id_ni) {
        var = pars->ncold * lnL->evaluatePartialAtP(pars->ir, 0, pars->id_unknown, pars->n); 
    } else if (pars->id_unknown == id_Tcold) {
        var = pars->ncold * lnL->evaluatePartialAtP(pars->ir, 0, pars->id_unknown, 0 /*n not used*/);
    } else {
        var = 0;
    }
    return prefactor * var;
}

/**
* Derivative of sp_cond (electric conductivity) wrt unknowns
*/
real_t HottailRateTermLowZ::dSpConddu(struct dGammadu_Params * pars){

    real_t res;
    if (pars->id_unknown == id_ni) {
        res = runawayFluid->evaluatePartialContributionBraamsConductivity(pars->ir, pars->id_unknown, pars->n);
    } else {
        res = 0;
    }
    return res;
}

/**
* Derivative of f0 wrt tau at p (integrand in other fcns)
*/
real_t HottailRateTermLowZ::df0dtau(real_t p, void *par){

    dGammadu_Params *pars = (dGammadu_Params*)par;
    
    real_t df0_dtau;
    pars->dist->evaluateEnergyDistributionFromTau(pars->ir, p, pars->tau, nullptr, nullptr, nullptr, &df0_dtau);
    
    return df0_dtau;
}

/**
* Integrand for the derivative (wrt tau) of the integral in the formula for Epar, needed in d/du (Ec/Epar^2) and other expressions
*/
real_t HottailRateTermLowZ::totalIntegrandFor_dBoxIntdu(real_t p, void *par){
    real_t totalIntegrand = df0dtau(p, par)*partialIntegrandForEpar(p);
    return totalIntegrand;
}

/**
* Derivative (wrt tau) of the integral in the formula for Epar, needed in d/du (Ec/Epar^2) and other expressions
*/
real_t HottailRateTermLowZ::dBoxIntdu(struct dGammadu_Params * pars){

    real_t result;
    if (pars->id_unknown == id_tau) { //f0 only depends on tau
    
        real_t error;
        gsl_function Func;
        Func.function = &totalIntegrandFor_dBoxIntdu;
        Func.params = pars;
        
        gsl_integration_qagiu(&Func, 0, ABSTOL_FOR_INT, RELTOL_FOR_INT, nGslIntervals, w, &result, &error);
        
    } else {
        result = 0;
    }
    
    return result;
}


/**
* Integral of df0/dtau from sqrt(Ec/Epar) to inf, part of the derivative of the integral in the expression for gamma
*/
real_t HottailRateTermLowZ::Intdf0du(struct dGammadu_Params * pars){

    real_t result;
    if (pars->id_unknown == id_tau) { //f0 only depends on tau
        
        real_t error;
        gsl_function Func;
        Func.function = &df0dtau;
        Func.params = pars;
        
        gsl_integration_qagiu(&Func, pars->lowlim, ABSTOL_FOR_INT, RELTOL_FOR_INT, nGslIntervals, w, &result, &error);
    
    } else {
        result = 0;
    }
    
    return result;
}



/**
* dgamma/du at ir (for ion and charge index n when id_unknown is id_ni)
*/
real_t HottailRateTermLowZ::dGammadu(len_t ir, len_t id_unknown, len_t n){

    // Create struct that holds relevant constants and values
    struct dGammadu_Params* derpar = new dGammadu_Params();
    
    derpar->ir = ir;
    derpar->id_unknown = id_unknown;
    derpar->ncold = unknowns->GetUnknownData(id_ncold)[ir];
    derpar->lnL = lnL->evaluateAtP(ir,0);
    derpar->tau = unknowns->GetUnknownData(id_tau)[ir];
    derpar->dist = distHT;
    derpar->n = n;
    
    real_t Ec = 4*M_PI* derpar->ncold * derpar->lnL *Constants::r0*Constants::r0*Constants::c*Constants::c*Constants::me / Constants::ec;
    real_t Zeff = ionHandler->GetZeff(ir);
    real_t sp_cond = runawayFluid->evaluateBraamsElectricConductivity(ir, T_final[ir], Zeff);
    real_t integPrefactor =  Constants::ec*Constants::c / (3*(1+Zeff));
    real_t j0 = unknowns->GetUnknownInitialData(id_johm)[ir];
    
    // Denominator in expression for Epar (integral only needs to be calculated once)
    real_t Box = j0*Ec / Epar[ir]; 
    
    // Deriv. of Ec
    real_t dEc_du = dEcdu(derpar);
    // Deriv. of sp_cond
    real_t dSpCond_du = dSpConddu(derpar);
    // Deriv. of f0 exists as integrand function
    
    // Derivative of 1/Ec wrt unknown u at radial point ir
    real_t dInvOfEc_du = - dEc_du / (Ec*Ec);   
    // Derivative of denominator ("Box") in expression for Epar wrt unknown u at radial point ir
    real_t dBox_du = dEc_du * sp_cond + dSpCond_du * Ec + integPrefactor * dBoxIntdu(derpar);
    
    // Deriv. of Ec/Epar^2 
    real_t dEcOverEparSq_du = Box*Box* dInvOfEc_du / (j0*j0) + 2*Box* dBox_du / (j0*j0*Ec);
    // Deriv. of Epar
    real_t dEpar_du = j0*dEc_du / Box - Ec*j0*dBox_du / (Box*Box);
    // Deriv. of dEpar/dt (since Epar_prev is constant)
    real_t ddEpardt_du = dEpar_du / dt;
    
    // Deriv of lower limit in gamma integral sqrt(Ec/Epar)
    real_t dlowlim_du = (dEc_du / sqrt(Ec*Epar[ir]) - sqrt(Ec)*dEpar_du / sqrt(Epar[ir]*Epar[ir]*Epar[ir]))/2;
    derpar->lowlim = sqrt(Ec/Epar[ir]);
    
    // Evaluate f0 at lowlim
    real_t f0_lowlim = derpar->dist->evaluateEnergyDistributionFromTau(ir,derpar->lowlim,derpar->tau);
    // Integral of df0/du from lowlim to inf
    real_t int_df0du = Intdf0du(derpar);
    // Deriv. of integral in gamma with u dependance in both f0 and lowlim, acc. to Leibniz integral rule
    real_t dGammaInt_du = -f0_lowlim * dlowlim_du + int_df0du;
    
    // dEpar/dt
    real_t dotEpar = (Epar[ir] - Epar_prev[ir]) / dt;
    if (dotEpar < 0){ // ensure non-negative runaway rate
        dotEpar = 0;
    }
    
    // Saved value of the integral in the expression for gamma
    real_t gammaIntegral = integraloff0inGamma[ir];
    // Deriv. of gamma wrt u, composed of sum from contributions by alpha, beta, delta
    real_t dGamma_du = dEcOverEparSq_du * dotEpar * gammaIntegral + ddEpardt_du * gammaIntegral * Ec / (Epar[ir]*Epar[ir]) + dGammaInt_du * dotEpar * Ec / (Epar[ir]*Epar[ir]);
    
    return dGamma_du;
    
}


/**
 * Sets the Jacobian of this equation term
 */
bool HottailRateTermLowZ::SetJacobianBlock(const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(!HasJacobianContribution(derivId))
        return false;

    for(len_t ir=0; ir<nr; ir++){
        real_t V = GetVolumeScaleFactor(ir);
        if(derivId == id_ni){
            len_t nZ  = ionHandler->GetNZ(); // Total number of ion types
            const len_t *Zs = ionHandler->GetZs(); // Highest charge state for each ion type in order of ion type index
            for(len_t iz = 0; iz<nZ; iz++){ // loop ion type index
                for(len_t Z0=0; Z0<=Zs[iz]; Z0++){ // loop charge states
                    len_t n = ionHandler->GetIndex(iz,Z0); // index in list [nD0 nD1 nNe0 ... nNe10] if D and Ne
                    real_t dGamma = dGammadu(ir, derivId, n);
                    jac->SetElement(ir, n*nr+ir, scaleFactor * dGamma * V);
                }
            }
        } else {
            real_t dGamma = dGammadu(ir, derivId, 0); // Independent of n
            jac->SetElement(ir, ir, scaleFactor * dGamma * V);
        }
    }
    return true;
}


/**
 * Deallocator
 */
void HottailRateTermLowZ::Deallocate(){
    if(Epar_prev != nullptr){
        delete [] Epar_prev;
        delete [] params;
        delete [] integraloff0inGamma;
    }
    if(Epar != nullptr){
        delete [] Epar;
    }
}



