/**
 * Description
 */

#include "DREAM/Equations/Kinetic/BoltzmannKnockOnOperator.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/DREAMException.hpp"

using namespace DREAM;

BoltzmannKnockOnOperator::BoltzmannKnockOnOperator(FVM::Grid *g, FVM::UnknownQuantityHandler *u) 
    : EquationTerm(g), unknowns(u), id_ntot(u->GetUnknownID(OptionConstants::UQTY_N_TOT))
{
    GridRebuilt();
}

BoltzmannKnockOnOperator::~BoltzmannKnockOnOperator(){
    Deallocate();
}

real_t BoltzmannKnockOnOperator::PitchIntegrandAtanFunc(
    real_t xi_lower, real_t xi_upper, real_t xiPrime, real_t xiStar,
    real_t theta_lower, real_t theta_upper, real_t thetaPrime, real_t thetaStar,
    pitch_interval interval
){
    real_t xi_u=0, xi_l=0;
    switch(interval){
        case PITCH_INTERVAL_LOWER: {
            xi_u = xi_upper;
            if(xiPrime >= cos(theta_lower + thetaStar) && xiPrime <= cos(theta_lower - thetaStar) )
                xi_l = xi_lower;
            else 
                xi_l = cos(thetaPrime + thetaStar);
        } break;
        case PITCH_INTERVAL_MIDDLE: {
            if(xiPrime <= cos(theta_lower - thetaStar))
                xi_l = xi_lower;
            else 
                xi_l = cos(thetaPrime + thetaStar);
            if(xiPrime >= cos(theta_upper + thetaStar))
                xi_u = xi_upper;
            else 
                xi_u = cos(thetaPrime - thetaStar);
        } break;
        case PITCH_INTERVAL_UPPER: {
            xi_l = xi_lower;
            if(xiPrime >= cos(theta_upper + thetaStar) && xiPrime <= cos(theta_upper - thetaStar))
                xi_u = xi_upper;
            else {
                xi_u = cos(thetaPrime - thetaStar);
            }
        } break;
        default:
            throw DREAMException("BoltzmannKnockOnOperator: Invalid poloidal angle interval; xi interval straddles xiStar.");
    }

    #define atanFunc(XI) atan( ((XI)-xiStar*xiPrime) \
        / sqrt( (1-xiStar*xiStar)*(1-xiPrime*xiPrime) \
            - ((XI)-xiStar*xiPrime)* ((XI)-xiStar*xiPrime) ) )

    return ( atanFunc(xi_u) - atanFunc(xi_l));
    #undef atanFunc
}

real_t BoltzmannKnockOnOperator::PitchIntegrandAtXi(
    real_t xi_lo, real_t xi_up, real_t xiPrime, real_t thetaPrime, real_t xiStar, real_t thetaStar
){
    real_t theta_lo = acos(xi_lo);
    real_t theta_up = acos(xi_up);

    real_t val=0;
    if(xi_lo < -xiStar && xi_up > -xiStar){
        val += PitchIntegrandAtanFunc(xi_lo, -xiStar, xiPrime, xiStar, theta_lo, M_PI-thetaStar, thetaPrime, thetaStar, PITCH_INTERVAL_LOWER);
        if(xi_up >= xiStar)
            val += PitchIntegrandAtanFunc(-xiStar, xiStar, xiPrime, xiStar, M_PI-thetaStar, thetaStar, thetaPrime, thetaStar, PITCH_INTERVAL_MIDDLE)
                 + PitchIntegrandAtanFunc(xiStar, xi_up, xiPrime, xiStar, thetaStar, theta_up, thetaPrime, thetaStar, PITCH_INTERVAL_UPPER);
        else 
            val += PitchIntegrandAtanFunc(-xiStar, xi_up, xiPrime, xiStar, M_PI-thetaStar, theta_up, thetaPrime, thetaStar, PITCH_INTERVAL_LOWER);
    } else if (xi_lo < xiStar && xi_up > xiStar){
        val += PitchIntegrandAtanFunc(xi_lo, xiStar, xiPrime, xiStar, thetaStar, theta_lo, thetaPrime,thetaStar,PITCH_INTERVAL_MIDDLE)
             + PitchIntegrandAtanFunc(xiStar, xi_up, xiPrime, xiStar,thetaStar,theta_up, thetaPrime,thetaStar,PITCH_INTERVAL_UPPER);
    } else {
        pitch_interval interval;
        if(xi_up <= -xiStar)
            interval = PITCH_INTERVAL_LOWER;
        else if(xi_up <= xiStar)
            interval = PITCH_INTERVAL_MIDDLE;
        else 
            interval = PITCH_INTERVAL_UPPER; 
        val += PitchIntegrandAtanFunc(xi_lo, xi_up, xiPrime, xiStar, theta_lo, theta_up, thetaPrime, thetaStar, interval);
    }
    return val;
}


real_t BoltzmannKnockOnOperator::EvaluatePitchIntegrandAtTheta(real_t theta, void *par){
    ParametersForPitchTerm *params = (ParametersForPitchTerm*) par;
    FVM::FluxSurfaceAverager *FSA = params->rGrid->GetFluxSurfaceAverager();

    len_t ir = params->ir;
    real_t B, Jacobian, ROverR0, NablaR2;
    FSA->GeometricQuantitiesAtTheta(ir, theta, B, Jacobian, ROverR0, NablaR2, FVM::FLUXGRIDTYPE_DISTRIBUTION);
    real_t Bmin = params->rGrid->GetBmin(ir);
    real_t Bmax = params->rGrid->GetBmax(ir);
    real_t BOverBmin = (Bmin != 0) ? B/Bmin : 1;
    real_t xi0Prime = params->xi0Prime;
    real_t xiPrime = xi0Prime*FVM::MomentumGrid::evaluateXiOverXi0(xi0Prime, BOverBmin);
    real_t thetaPrime = acos(xiPrime);
    real_t preFactor = Jacobian*BOverBmin*xi0Prime/xiPrime;

    real_t xiStar = params->xiStar;
    real_t thetaStar = params->thetaStar;

    // TODO: identify trapped region and add mirrored intervals
    real_t xi0_up = params->xi0_jp;
    real_t xi0_lo = params->xi0_jm;
    real_t xi0Trapped = params->rGrid->GetXi0TrappedBoundary(ir);
    
    // entirely passing
    if(xi0_lo >= xi0Trapped || xi0_up <= -xi0Trapped){
        real_t xi_up = xi0_up*FVM::MomentumGrid::evaluateXiOverXi0(xi0_up, BOverBmin);
        real_t xi_lo = xi0_lo*FVM::MomentumGrid::evaluateXiOverXi0(xi0_lo, BOverBmin);
        return preFactor*PitchIntegrandAtXi(xi_lo, xi_up, xiPrime, thetaPrime, xiStar, thetaStar);
    } // otherwise, continue and handle the trapped case
    
    real_t val = 0;
    
    // local pitch of a particle on the (positive) trapped-passing boundary
    real_t xiTAtTheta = sqrt(1-B/Bmax); 

    // Pitch interval straddles negative trapped-passing boundary:
    // add the contribution from the passing-region part
    if(xi0_lo < -xi0Trapped){
        real_t xi_up = -xiTAtTheta;
        real_t xi_lo = xi0_lo*FVM::MomentumGrid::evaluateXiOverXi0(xi0_lo, BOverBmin);
        val += preFactor*PitchIntegrandAtXi(xi_lo, xi_up, xiPrime, thetaPrime, xiStar, thetaStar);
        if(xi0_up <=0) // no contribution from negative trapped region
            return val;
        xi0_lo = 0;
    } else if(xi0_lo < 0){
        if(xi0_up <= 0) // no contribution from negative trapped region
            return 0;
        else 
            xi0_lo = 0;
    }
    // Pitch interval straddles positive trapped-passing boundary:
    // add the contribution from the passing-region part
    if(xi0_up > xi0Trapped){
        real_t xi_up = xi0_up*FVM::MomentumGrid::evaluateXiOverXi0(xi0_up, BOverBmin);
        real_t xi_lo = xiTAtTheta;
        val += preFactor*PitchIntegrandAtXi(xi_lo, xi_up, xiPrime, thetaPrime, xiStar, thetaStar);
        xi0_up = xi0Trapped;        
    }

    xi0_lo = std::max( xiTAtTheta, xi0_lo );
    real_t xi_lo = xi0_lo*FVM::MomentumGrid::evaluateXiOverXi0(xi0_lo, BOverBmin);
    real_t xi_up = xi0_up*FVM::MomentumGrid::evaluateXiOverXi0(xi0_up, BOverBmin);
    
    // contribution from trapped orbit including its mirror
    val += preFactor * (
        PitchIntegrandAtXi(xi_lo, xi_up, xiPrime, thetaPrime, xiStar, thetaStar)
        + PitchIntegrandAtXi(-xi_up, -xi_lo, xiPrime, thetaPrime, xiStar, thetaStar)
    );

    return val;
}


void BoltzmannKnockOnOperator::SetVectorElements(real_t *vec, const real_t *f){
    const real_t *ntot = unknowns->GetUnknownData(id_ntot);
    for(len_t ir=0; ir<nr; ir++){
        len_t Nind = n1[ir]*n2[ir];
        for(len_t i=0; i<Nind; i++){
            len_t ind = ir*Nind+i;
            for(len_t j=0; j<Nind; j++)
                vec[ind] += ntot[ir]*SourceMatrix[ir][i][j]*f[ir*Nind+j];
        }
    }   
}

bool BoltzmannKnockOnOperator::GridRebuilt(){
    Deallocate();
    this->EquationTerm::GridRebuilt();
    Allocate();

    return true;
}

/**
 * Evaluates the energy-dependent function appearing in the Möller cross-section, 
 * corresponding to the function Sigma given by Eq (5) in doc/notes/knockon
 *   gamma: Lorentz factor of the knock-on electron
 *  gamma1: Lorentz factor of the incident primary electron
 */
real_t BoltzmannKnockOnOperator::EvaluateMollerDifferentialCrossSection(real_t gamma, real_t gamma1) {
    real_t gamma1Sq = gamma1*gamma1;
    real_t s = (gamma-1)*(gamma1-gamma);
    real_t preFactor = 2*M_PI*Constants::r0 * Constants::r0 * Constants::c
        * gamma1Sq / ( (gamma1Sq-1) * s * s );
    return preFactor * (
        (gamma1-1)*(gamma1-1) - s/gamma1Sq * (
            2*gamma1Sq + 2*gamma1 - 1 - s
        )
    );
}

/**
 * Evaluates v1 = p1/gamma1 times the integral over gamma of 'EvaluateMollerDifferentialCrossSection'
 * between gamma=gamma_l and gamma_u. 
 * Corresponds to the Sigma_i function given by Eq (11) in doc/notes/knockon.
 */
real_t BoltzmannKnockOnOperator::EvaluateIntegratedMollerDCS(real_t gamma_l, real_t gamma_u, real_t gamma1){
    real_t 
        p1 = sqrt(gamma1*gamma1-1),
        preFactor = 2*M_PI*Constants::r0*Constants::r0*Constants::c/(gamma1*p1),
        lnTerm_u = log( (gamma1-gamma_u)/(gamma_u-1) ),
        lnTerm_l = log( (gamma1-gamma_l)/(gamma_l-1) ),
        g1Sq = gamma1*gamma1,
        pre1 = g1Sq,
        pre2 = -(g1Sq + 2*gamma1 - 1),
        pre3 = 1.0;
    
    real_t Term1 = 1/(gamma1-gamma_u) - 1/(gamma1-gamma_l) - 1/(gamma_u-1) + 1/(gamma_l-1);
    real_t Term2 = - 1/(gamma1-1) * (lnTerm_u - lnTerm_l);
    real_t Term3 = gamma_u - gamma_l;

    return preFactor * ( pre1*Term1 + pre2*Term2 + pre3*Term3 );
}

/**
 * Essentially full rows in momentum are possible. In practice most 
 * rows will have significantly less than half of these elements
 * assigned, but it varies dramatically from point to point
 */
len_t BoltzmannKnockOnOperator::GetNumberOfNonZerosPerRow() const {
    return grid->GetMaxNp1()*grid->GetMaxNp2();
}
