#include "DREAM/Equations/Fluid/ExternalAvalancheTerm.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"

/**
 * Implementation of an equation term that matches the 
 * 'fluid' mode of AvalancheSourceRP, which integrates
 * the RP source over some p>pCutoff, to the analytic
 * growth rate of AvalancheGrowthTerm, using the formula
 *      G = (G1^a + G2^a)^(1/a),
 * where for a<0 is closest to the smallest of the two.
 * 
 * The purpose is to resolve the issue where the external
 * avalanche rate of AvalancheSourceRP is non-zero positive
 * for all electric fields, which crashes the solution
 * when you enter a situation where E<Eceff and you should 
 * have decay. In this case, this equation term will be 
 * exactly the AvalancheGrowthTerm.
 */

using namespace DREAM;

/**
 * Constructor
 */
ExternalAvalancheTerm::ExternalAvalancheTerm(
    FVM::Grid *g, real_t pc, real_t a, RunawayFluid *ref, 
    FVM::UnknownQuantityHandler *u, real_t sf
) : FVM::DiagonalComplexTerm(g,u), REFluid(ref), pCutoff(pc), aExp(a), scaleFactor(sf) {
    id_ntot   = unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT);
    id_Efield = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    
    AddUnknownForJacobian(unknowns, id_Efield);
    AddUnknownForJacobian(unknowns, id_ntot);
}

/**
 * Sets the weights of this diagonal equation term
 */ 
void ExternalAvalancheTerm::SetWeights(){
    const real_t *GammaFluid = REFluid->GetAvalancheGrowthRate();
    const real_t *n_tot = unknowns->GetUnknownData(id_ntot);

    for(len_t ir=0; ir<nr; ir++){
        real_t GammaExternal = n_tot[ir]*AvalancheSourceRP::EvaluateNormalizedTotalKnockOnNumber(
            grid->GetRadialGrid()->GetFSA_B(ir), pCutoff
        );
        real_t GammaMatched = GammaFluid[ir];
        if(GammaMatched>0)
            GammaMatched = pow( pow(GammaExternal,aExp)  + pow(GammaFluid[ir],aExp) , 1.0/aExp);
        weights[ir] = scaleFactor * GammaMatched;
    }
}

/**
 * Sets the jacobian of the weights of this equation term.
 */
void ExternalAvalancheTerm::SetDiffWeights(len_t derivId, len_t /*nMultiples*/) {
    AllocateDGamma();
    REFluid->evaluatePartialContributionAvalancheGrowthRate(dGammaFluid, derivId);
    const real_t *GammaFluid = REFluid->GetAvalancheGrowthRate();
    const real_t *n_tot = unknowns->GetUnknownData(id_ntot);
    for(len_t ir=0; ir<nr; ir++){
        if(GammaFluid[ir]<=0)
            diffWeights[ir] = scaleFactor * dGammaFluid[ir];
        else {
            real_t GammaExternal = n_tot[ir]*AvalancheSourceRP::EvaluateNormalizedTotalKnockOnNumber(
                grid->GetRadialGrid()->GetFSA_B(ir), pCutoff
            );
            real_t dGammaExternal = 0;
            if(derivId==id_ntot) // add external-contribution jacobian
                dGammaExternal = AvalancheSourceRP::EvaluateNormalizedTotalKnockOnNumber(
                    grid->GetRadialGrid()->GetFSA_B(ir), pCutoff
                );
            real_t Factor1 = pow(GammaExternal,aExp-1)*dGammaExternal + pow(GammaFluid[ir],aExp-1)*dGammaFluid[ir]; 
            real_t Factor2 = pow(  pow(GammaExternal,aExp) + pow(GammaFluid[ir],aExp) , 1.0/aExp - 1.0);
            diffWeights[ir] = scaleFactor * Factor1*Factor2;
        }
    }
}
