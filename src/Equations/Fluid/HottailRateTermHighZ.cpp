/**
 * Implementation of equation term representing the runaway generation rate
 * due to hottail when using an analytic distribution function
 */

#include "DREAM/Equations/Fluid/HottailRateTermHighZ.hpp"

using namespace DREAM;

/**
 * Constructor.
 * 
 * gsl_altPc* contains parameters and functions needed
 * to evaluate the critical runaway momentum in the hottail
 * calculation using a gsl root-finding algorithm
 * (using the 'alternative' model for pc in Ida's MSc thesis)
 */
HottailRateTermHighZ::HottailRateTermHighZ(           
    FVM::Grid *grid, AnalyticDistributionHottail *dist, FVM::UnknownQuantityHandler *unknowns,
    IonHandler *ionHandler, CoulombLogarithm *lnL, real_t sf
) : HottailRateTerm(grid, dist, unknowns,sf), lnL(lnL),
    id_ncold(unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD)),
    id_Efield(unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD)),
    id_tau(unknowns->GetUnknownID(OptionConstants::UQTY_TAU_COLL))
{
    SetName("HottailRateTermHighZ");

    AddUnknownForJacobian(unknowns,id_Efield);
    AddUnknownForJacobian(unknowns,id_ncold);
    AddUnknownForJacobian(unknowns,id_tau);
    //AddUnknownForJacobian(unknowns,id_ni);    // Zeff and lnL (nfree) jacobian
    //AddUnknownForJacobian(unknowns,id_Tcold); // lnL jacobian
    

    this->fdfsolver = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_secant);

    gsl_params.ionHandler = ionHandler;
    gsl_params.rGrid = grid->GetRadialGrid();
    gsl_params.dist = dist;

    gsl_func.f = &(PcFunc);
    gsl_func.df = &(PcFunc_df);
    gsl_func.fdf = &(PcFunc_fdf);
    gsl_func.params = &gsl_params;

    this->GridRebuilt();
}

/**
 * Destructor
 */
HottailRateTermHighZ::~HottailRateTermHighZ(){
    Deallocate();
    gsl_root_fdfsolver_free(fdfsolver);
}           

/**
 * Called after the grid is rebuilt; (re)allocates memory
 * for all quantities
 */
bool HottailRateTermHighZ::GridRebuilt(){
    this->HottailRateTerm::GridRebuilt();

    Deallocate();
    pCrit_prev = new real_t[nr];

    return true;
}

/**
 * Rebuilds quantities used by this equation term.
 * Note that the equation term uses the _energy distribution_
 * from AnalyticDistributionHottail which differs from the 
 * full distribution function by a factor of 4*pi
 */
void HottailRateTermHighZ::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler*) {
    this->dt = dt;
    bool newTimeStep = (t!=tPrev);
    if(newTimeStep)
        tPrev = t;
    for(len_t ir=0; ir<nr; ir++){
        if(newTimeStep)
            pCrit_prev[ir] = pCrit[ir];

        real_t fAtPc, dfdpAtPc;
        pCrit[ir] = evaluateCriticalMomentum(ir, fAtPc, dfdpAtPc);
        real_t dotPc = (pCrit[ir] - pCrit_prev[ir]) / dt;
        if (dotPc > 0) // ensure non-negative runaway rate
            dotPc = 0;
        gamma[ir] = -pCrit[ir]*pCrit[ir]*dotPc*fAtPc; // generation rate
        // set derivative of gamma with respect to pCrit (used for jacobian)
        dGammaDPc[ir] = -(2*pCrit[ir]*dotPc*fAtPc + pCrit[ir]*pCrit[ir]*fAtPc/dt + pCrit[ir]*pCrit[ir]*dotPc*dfdpAtPc);
    }
}

/**
 * Function whose root (with respect to p) represents the 
 * critical runaway momentum in the 'alternative' model
 */
real_t HottailRateTermHighZ::PcFunc(real_t p, void *par) { 
    if(p<0) // handle case where algorithm leaves the physical domain of non-negative momenta
        p=0;
    PcParams *params = (PcParams*)par;

    len_t ir = params->ir;
    real_t Eterm = params->Eterm;
    real_t ncold = params->ncold;
    real_t tau   = params->tau;
    real_t lnL = params->lnL;
    real_t dFdpOverF;
    params->F = params->dist->evaluateEnergyDistributionFromTau(ir,p,tau,&params->dFdp,nullptr, &dFdpOverF);

    real_t Ec = 4*M_PI*ncold*lnL*Constants::r0*Constants::r0*Constants::c * Constants::me * Constants::c / Constants::ec;
    real_t E = Eterm/Ec;
    real_t EPF = params->rGrid->GetEffPassFrac(ir);
    real_t Zeff = params->ionHandler->GetZeff(ir);

    real_t p2 = p*p;
    real_t gamma = sqrt(1+p2); 
     
    // for the non-relativistic distribution, this function is
    // approximately linear, yielding efficient root finding
    return  sqrt((p/gamma)*cbrt( p2*E*E*EPF * (-dFdpOverF) )) - sqrt(cbrt( 3*(1+Zeff)));
    // previous equivalent expression:
    // real_t g3 = (1+p2)*gamma;
    // return  sqrt(cbrt( p2*p2*p*E*E*EPF * (-dFdpOverF) )) - sqrt(cbrt( 3.0*(1+Zeff)*g3));
}

/**
 * Returns the derivative of PcFunc with respect to p
 */
real_t HottailRateTermHighZ::PcFunc_df(real_t p, void *par) {
    real_t h = 1e-3*p;
    return (PcFunc(p+h,par) - PcFunc(p,par)) / h;
}

/**
 * Method which sets both f=PcFunc and df=PcFunc_df
 */
void HottailRateTermHighZ::PcFunc_fdf(real_t p, void *par, real_t *f, real_t *df){
    real_t h = 1e-3*p;
    *f = PcFunc(p,par);
    *df = (PcFunc(p+h, par) - *f) / h;
}

/**
 * Evaluates the 'alternative' critical momentum pc using Ida's MSc thesis (4.35) 
 */
real_t HottailRateTermHighZ::evaluateCriticalMomentum(len_t ir, real_t &f, real_t &dfdp){
    gsl_params.ir = ir;
    gsl_params.lnL   = lnL->evaluateAtP(ir,0);
    gsl_params.ncold = unknowns->GetUnknownData(id_ncold)[ir];
    gsl_params.Eterm = unknowns->GetUnknownData(id_Efield)[ir];
    gsl_params.tau   = unknowns->GetUnknownData(id_tau)[ir];
    
    real_t root = (pCrit_prev[ir] == 0) ? 5*distHT->GetInitialThermalMomentum(ir) : pCrit_prev[ir];
    RunawayFluid::FindRoot_fdf_bounded(0,std::numeric_limits<real_t>::infinity(),root, gsl_func, fdfsolver, RELTOL_FOR_PC, ABSTOL_FOR_PC);
    f = gsl_params.F;
    dfdp = gsl_params.dFdp;
    return root;
}

/**
 * Evaluates the jacobian of CriticalMomentum with 
 * respect to the unknown with id 'derivId'
 */
real_t HottailRateTermHighZ::evaluatePartialCriticalMomentum(len_t ir, len_t derivId){
    gsl_params.ir = ir;
    gsl_params.lnL   = lnL->evaluateAtP(ir,0);
    gsl_params.ncold = unknowns->GetUnknownData(id_ncold)[ir];
    gsl_params.Eterm = unknowns->GetUnknownData(id_Efield)[ir];
    gsl_params.tau   = unknowns->GetUnknownData(id_tau)[ir];

    real_t h = 0;
    if(derivId == id_Efield){
        h = (1.0+fabs(gsl_params.Eterm))*sqrt(RELTOL_FOR_PC),
        gsl_params.Eterm += h;
    } else if (derivId == id_ncold){
        h = (1.0+fabs(gsl_params.ncold))*sqrt(RELTOL_FOR_PC),
        gsl_params.ncold += h;
    } else if (derivId == id_tau){
        h = (1.0+fabs(gsl_params.tau))*sqrt(RELTOL_FOR_PC),
        gsl_params.tau += h;
    }

    real_t root = pCrit[ir];
    RunawayFluid::FindRoot_fdf_bounded(0,std::numeric_limits<real_t>::infinity(),root, gsl_func, fdfsolver, RELTOL_FOR_PC, ABSTOL_FOR_PC);
    return (root - pCrit[ir])/h;
}


/**
 * Sets the Jacobian of this equation term
 */
bool HottailRateTermHighZ::SetJacobianBlock(const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(!HasJacobianContribution(derivId))
        return false;

    for(len_t ir=0; ir<nr; ir++){
        const len_t xiIndex = this->GetXiIndexForEDirection(ir);
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
        real_t V = GetVolumeScaleFactor(ir);

        // Check if the quantity w.r.t. which we differentiate is a
        // fluid quantity, in which case it has np1=1, xiIndex=0
        len_t np1_op = np1, xiIndex_op = xiIndex;
        if (unknowns->GetUnknown(derivId)->NumberOfElements() == nr) {
            np1_op = 1;
            xiIndex_op = 0;
        }
        real_t dPc = evaluatePartialCriticalMomentum(ir, derivId);
        real_t dGamma = dPc * dGammaDPc[ir];

        if(derivId==id_tau){ // add contribution from explicit tau dependence in f
            real_t dotPc = (pCrit[ir] - pCrit_prev[ir]) / dt;
            if (dotPc > 0) // ensure non-negative runaway rate
                dotPc = 0;
            real_t tau = unknowns->GetUnknownData(id_tau)[ir];
            real_t dFdTauAtPc;
            distHT->evaluateEnergyDistributionFromTau(ir, pCrit[ir], tau, nullptr, nullptr, nullptr, &dFdTauAtPc);
            dGamma -= pCrit[ir]*pCrit[ir]*dotPc*dFdTauAtPc;
        }
        jac->SetElement(ir + np1*xiIndex, ir + np1_op*xiIndex_op, scaleFactor * dGamma * V);
    }

    return true;
}

/**
 * Deallocator
 */
void HottailRateTermHighZ::Deallocate(){
    if(pCrit_prev != nullptr)
        delete [] pCrit_prev;
}
