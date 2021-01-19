/**
 * Implementation of equation term representing the runaway generation rate
 * due to hottail when using an analytic distribution function
 */

#include "DREAM/Equations/Fluid/HottailRateTerm.hpp"

using namespace DREAM;

/**
 * Constructor.
 * 
 * gsl_altPc* contains parameters and functions needed
 * to evaluate the critical runaway momentum in the hottail
 * calculation using a gsl root-finding algorithm
 * (using the 'alternative' model for pc in Ida's MSc thesis)
 */
HottailRateTerm::HottailRateTerm(           
    FVM::Grid *grid, AnalyticDistribution *dist, FVM::UnknownQuantityHandler *unknowns,
    IonHandler *ionHandler, CoulombLogarithm *lnL,
    enum OptionConstants::eqterm_hottail_mode type, real_t sf
) : FVM::EquationTerm(grid), RunawaySourceTerm(grid, unknowns), type(type), scaleFactor(sf), distHT(dist), unknowns(unknowns), lnL(lnL),
    id_ncold(unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD)),
    id_Efield(unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD))
{
    this->fdfsolver = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_secant);

    gsl_altPcParams.ionHandler = ionHandler;
    gsl_altPcParams.rGrid = grid->GetRadialGrid();
    gsl_altPcParams.dist = dist;

    gsl_altPcFunc.f = &(altPcFunc);
    gsl_altPcFunc.df = &(altPcFunc_df);
    gsl_altPcFunc.fdf = &(altPcFunc_fdf);
    gsl_altPcFunc.params = &gsl_altPcParams;

    this->GridRebuilt();
}

/**
 * Destructor
 */
HottailRateTerm::~HottailRateTerm(){
    DeallocateAll();
    gsl_root_fdfsolver_free(fdfsolver);
}           

/**
 * Called after the grid is rebuilt; (re)allocates memory
 * for all quantities
 */
bool HottailRateTerm::GridRebuilt(){
    this->EquationTerm::GridRebuilt();

    DeallocateAll();
    pcAlt      = new real_t[nr];
    pcAlt_prev = new real_t[nr];
    gamma      = new real_t[nr];
    dGammaDPc  = new real_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        pcAlt[ir] = 0;
    }

    return true;
}

/**
 * Rebuilds quantities used by this equation term
 */
void HottailRateTerm::Rebuild(const real_t, const real_t dt, FVM::UnknownQuantityHandler*) {
    this->dt = dt;
    if(type == OptionConstants::EQTERM_HOTTAIL_MODE_ANALYTIC_ALT_PC){ // Ida MSc thesis (4.39)
        for(len_t ir=0; ir<nr; ir++){
            pcAlt_prev[ir] = pcAlt[ir];
            pcAlt[ir] = evaluateAltCriticalMomentum(ir);
            real_t dfdp;
            real_t f = distHT->evaluateEnergyDistribution(ir,pcAlt[ir], &dfdp);
            real_t dotPcAlt = (pcAlt[ir] - pcAlt_prev[ir]) / dt;
            gamma[ir] = -4*M_PI*pcAlt[ir]*pcAlt[ir]*dotPcAlt*f; // generation rate
            // set derivative of gamma with respect to pcAlt (used for jacobian)
            dGammaDPc[ir] = -4*M_PI*(2*pcAlt[ir]*dotPcAlt*f + pcAlt[ir]*pcAlt[ir]*f/dt + pcAlt[ir]*pcAlt[ir]*dotPcAlt*dfdp);
        }
    } else if (type == OptionConstants::EQTERM_HOTTAIL_MODE_ANALYTIC) { // Ida MSc Thesis (4.24)
        // TODO
    }
}

/**
 * Function whose root (with respect to p) represents the 
 * critical runaway momentum in the 'alternative' model
 */
real_t HottailRateTerm::altPcFunc(real_t p, void *par) { 
    altPcParams *params = (altPcParams*)par;

    len_t ir = params->ir;
    real_t Eterm = params->Eterm;
    real_t ncold = params->ncold;
    real_t lnL = params->lnL;
    real_t dFdp;
    real_t F = params->dist->evaluateEnergyDistribution(ir,p,&dFdp);

    real_t Ec = 4*M_PI*ncold*lnL*Constants::r0*Constants::r0*Constants::c;
    real_t E = Eterm/Ec;
    real_t EPF = params->rGrid->GetEffPassFrac(ir);
    real_t Zeff = params->ionHandler->GetZeff(ir);

    real_t p2 = p*p;
//    return p2*p2*p* E*E*EPF * dFdp/F + 3.0*(1+Zeff);
    return sqrt(sqrt( p2*p2*p*E*E*EPF * (-dFdp/F) )) - sqrt(sqrt( 3.0*(1+Zeff) ));
}

/**
 * Returns the derivative of altPcFunc with respect to p
 */
real_t HottailRateTerm::altPcFunc_df(real_t p, void *par) {
    real_t h = 0.01*p;
    return (altPcFunc(p+h,par) - altPcFunc(p,par)) / h;
}

/**
 * Method which sets both f=altPcFunc and df=altPcFunc_df
 */
void HottailRateTerm::altPcFunc_fdf(real_t p, void *par, real_t *f, real_t *df){
    real_t h = 0.01*p;
    *f = altPcFunc(p,par);
    *df = (altPcFunc(p+h, par) - *f) / h;
}

/**
 * Evaluates the 'alternative' critical momentum pc using Ida's MSc thesis (4.35) 
 */
real_t HottailRateTerm::evaluateAltCriticalMomentum(len_t ir){
    real_t root = (pcAlt_prev[ir] == 0) ? 1 : pcAlt_prev[ir];
    gsl_altPcParams.ir = ir;
    gsl_altPcParams.lnL = lnL->evaluateAtP(ir,0);
    gsl_altPcParams.ncold = unknowns->GetUnknownData(id_ncold)[ir];
    gsl_altPcParams.Eterm = unknowns->GetUnknownData(id_Efield)[ir];

    RunawayFluid::FindRoot_fdf(root, gsl_altPcFunc, fdfsolver);
    return root;
}

/**
 * Evaluates the jacobian of AltCriticalMomentum with 
 * respect to the unknown with id 'derivId'
 */
real_t HottailRateTerm::evaluatePartialAltCriticalMomentum(len_t ir, len_t derivId){
    real_t root = pcAlt[ir];
    gsl_altPcParams.ir = ir;
    gsl_altPcParams.lnL = lnL->evaluateAtP(ir,0);
    gsl_altPcParams.ncold = unknowns->GetUnknownData(id_ncold)[ir];
    gsl_altPcParams.Eterm = unknowns->GetUnknownData(id_Efield)[ir];

    real_t h = 0;
    if(derivId == id_Efield){
        h = gsl_altPcParams.Eterm*.01,
        gsl_altPcParams.Eterm += h;
    } else if (derivId == id_ncold){
        h = gsl_altPcParams.ncold*.01,
        gsl_altPcParams.ncold += h;
    }

   RunawayFluid::FindRoot_fdf(root, gsl_altPcFunc, fdfsolver);
   return (root - pcAlt[ir])/h;
}

/**
 * Sets the matrix elements of this equation term
 */
void HottailRateTerm::SetMatrixElements(FVM::Matrix */*mat*/, real_t *rhs){
    SetVectorElements(rhs, nullptr);
}

/**
 * Sets the vector elements of this equation term
 */
void HottailRateTerm::SetVectorElements(real_t *vec, const real_t *){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        const len_t xiIndex = this->GetXiIndexForEDirection(ir);
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
        const len_t np2 = this->grid->GetMomentumGrid(ir)->GetNp2();
        real_t V = GetVolumeScaleFactor(ir);

        // Insert at p=0, xi=+1 (or -1 if E is negative)
        vec[offset + np1*xiIndex + 0] += scaleFactor*gamma[ir]*V;

        offset += np1*np2;
    }
}

/**
 * Sets the Jacobian of this equation term
 */
void HottailRateTerm::SetJacobianBlock(const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(!HasJacobianContribution(derivId))
        return;

    for(len_t ir=0; ir<nr; ir++){
        const len_t xiIndex = this->GetXiIndexForEDirection(ir);
        const len_t np1 = this->grid->GetMomentumGrid(ir)->GetNp1();
        real_t V = GetVolumeScaleFactor(ir);

        // Check if the quantity w.r.t. which we differentiate is a
        // fluid quantity, in which case it has np1=0, xiIndex=0
        len_t np1_op = np1, xiIndex_op = xiIndex;
        if (unknowns->GetUnknown(derivId)->NumberOfElements() == nr) {
            np1_op = 1;
            xiIndex_op = 0;
        }
        real_t dGamma;
        if(type == OptionConstants::EQTERM_HOTTAIL_MODE_ANALYTIC_ALT_PC){
            real_t dPc = evaluatePartialAltCriticalMomentum(ir, derivId);
            dGamma = dPc * dGammaDPc[ir];
        } else if ( type == OptionConstants::EQTERM_HOTTAIL_MODE_ANALYTIC) {
            // TODO
        }
        jac->SetElement(ir + np1*xiIndex, ir + np1_op*xiIndex_op, scaleFactor * dGamma * V);
    }
}

/**
 * Deallocator
 */
void HottailRateTerm::DeallocateAll(){
    if(pcAlt != nullptr){
        delete [] pcAlt;
        delete [] pcAlt_prev;
        delete [] gamma;
        delete [] dGammaDPc;
    }
}
