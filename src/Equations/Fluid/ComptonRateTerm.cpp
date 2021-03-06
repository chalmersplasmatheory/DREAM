#include "DREAM/Equations/Fluid/ComptonRateTerm.hpp"

/**
 * Implementation of a class which represents the Compton runaway rate contribution to the n_re equation.
 * Employs the analytical generation rate calculated by RunawayFluid.
 */
using namespace DREAM;

/**
 * Constructor.
 */
ComptonRateTerm::ComptonRateTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *uqn,
    RunawayFluid *rf, FVM::Grid *operandGrid, real_t scaleFactor
) : FVM::DiagonalComplexTerm(g,uqn, operandGrid), RunawaySourceTerm(g,uqn),
    REFluid(rf), scaleFactor(scaleFactor) {

    SetName("ComptonRateTerm");

    AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD));
    AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT));
}


/**
 * Destructor.
 */
ComptonRateTerm::~ComptonRateTerm() {
    if(dGamma != nullptr)
        delete [] dGamma;
}


/**
 * if nr has changed, (re)allocates memory for jacobian weights
 */
void ComptonRateTerm::AllocateDGamma(){
    if(nr_tmp != this->nr){
        if(dGamma != nullptr)
            delete [] dGamma;
        dGamma = new real_t[this->nr];
        nr_tmp = this->nr;
    }
}


/**
 * Sets jacobian of the weights of this diagonal term to the
 * derivative of the Compton runaway generation rate with respect
 * to the unknown with id derivId.
 */
void ComptonRateTerm::SetDiffWeights(len_t derivId, len_t nMultiples){
    AllocateDGamma();
    REFluid->evaluatePartialContributionComptonGrowthRate(dGamma, derivId);
    len_t offset = 0;
    for(len_t n = 0; n<nMultiples; n++)
        for (len_t ir = 0; ir < nr; ir++){
            const len_t xiIndex = this->GetXiIndexForEDirection(ir);
            const real_t V = this->GetVolumeScaleFactor(ir);

            diffWeights[offset + n1[ir]*xiIndex + 0] = scaleFactor*dGamma[ir] * V;
            offset += n1[ir]*n2[ir];
        }
}


/**
 * Sets weights of this diagonal term to the Compton runaway generation rate
 */
void ComptonRateTerm::SetWeights(){
    const real_t *comptonRate = REFluid->GetComptonRunawayRate();
    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++){
        const len_t xiIndex = this->GetXiIndexForEDirection(ir);
        const real_t V = this->GetVolumeScaleFactor(ir);

        weights[offset + n1[ir]*xiIndex + 0] = scaleFactor*comptonRate[ir] * V;
        offset += n1[ir]*n2[ir];
    }
}
