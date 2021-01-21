/**
 * Implementation of class representing the analytic f_hot distribution
 */

#include "DREAM/Equations/AnalyticDistributionHottail.hpp"
#include "DREAM/DREAMException.hpp"

using namespace DREAM;

/**
 * Constructor
 */
AnalyticDistributionHottail::AnalyticDistributionHottail(FVM::RadialGrid *rGrid, FVM::UnknownQuantityHandler *u, real_t *n0, real_t *T0, OptionConstants::uqty_f_hot_dist_mode type)
    : AnalyticDistribution(rGrid,u), type(type), n0(n0), T0(T0)
{
    if(u->HasUnknown(OptionConstants::UQTY_TAU_COLL))
        id_tau = u->GetUnknownID(OptionConstants::UQTY_TAU_COLL);

    GridRebuilt();
}

/**
 * Destructor
 */
AnalyticDistributionHottail::~AnalyticDistributionHottail(){
    if(preFactor != nullptr){
        delete [] preFactor;
        delete [] betaTh;
    }
    if(n0 !=nullptr){
        delete [] n0;
        delete [] T0;
    }
}

/**
 * Called when grid has been rebuilt: (re)allocates memory for all quantities 
 * and carries out calculations of constants
 */
bool AnalyticDistributionHottail::GridRebuilt(){
    this->AnalyticDistribution::GridRebuilt();

    if(preFactor != nullptr){
        delete [] preFactor;
        delete [] betaTh;
    }

    preFactor = new real_t[nr];
    betaTh    = new real_t[nr];
    for(len_t ir=0; ir<nr; ir++){
        betaTh[ir] = sqrt(2*T0[ir] / Constants::mc2inEV);
        preFactor[ir] = n0[ir]/(M_PI*M_SQRTPI*betaTh[ir]*betaTh[ir]*betaTh[ir]);
    }

    return true;
}

/**
 * Evaluate the pitch-angle averaged energy distribution
 */
real_t AnalyticDistributionHottail::evaluateEnergyDistribution(len_t ir, real_t p, real_t *dfdp, real_t *dfdr){
    if(type == OptionConstants::UQTY_F_HOT_DIST_MODE_NONREL){
        real_t tau = unknowns->GetUnknownData(id_tau)[ir];
        return evaluateEnergyDistributionFromTau(ir,p,tau,dfdp,dfdr);
    } else {
        return NAN;
        throw DREAMException("AnalyticDistributionHottail: Invalid type %d", type);
    }
}

/**
 * Evaluation of the analytic energy distribution given 
 * in Smith & Verwichte (2008), equations (9-10),
 * using the time-integrated (ideal) slowing-down frequency 'tau'
 */
real_t AnalyticDistributionHottail::evaluateEnergyDistributionFromTau(len_t ir, real_t p, real_t tau, real_t *dfdp, real_t *dfdr, real_t *dFdpOverF, real_t *dFdTau){
    if(type != OptionConstants::UQTY_F_HOT_DIST_MODE_NONREL)
        throw DREAMException("AnalyticDistributionHottail: Invalid type %d", type);

    real_t s = std::pow(p*p*p + 3*tau,2.0/3.0);
    real_t exponent = - s / (betaTh[ir]*betaTh[ir]);
    real_t f = preFactor[ir] * exp(exponent);

    if(dfdr != nullptr){
        if(nr==1)
            *dfdr = 0;
        else {
            len_t ind_shift = (ir<nr-1) ? ir+1 : ir-1; 
            real_t f_shift = evaluateEnergyDistributionFromTau(ind_shift, p, tau);
            real_t dr = rGrid->GetR(ir) - rGrid->GetR(ind_shift);
            *dfdr = (f-f_shift) / dr;
        }
    }

    if(dfdp != nullptr)
        *dfdp = - f * 2*p*p / ( sqrt(s)*betaTh[ir]*betaTh[ir] );

    if(dFdpOverF != nullptr)
        *dFdpOverF = - 2*p*p / ( sqrt(s)*betaTh[ir]*betaTh[ir] );

    if(dFdTau != nullptr)
        *dFdTau = - 2*f / ( sqrt(s)*betaTh[ir]*betaTh[ir] );

    return f;
}
