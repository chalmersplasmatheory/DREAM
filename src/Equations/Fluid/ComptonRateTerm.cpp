#include "DREAM/Equations/Fluid/ComptonRateTerm.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "DREAM/DREAMException.hpp"

/**
 * Implementation of a class which represents the Gamma_compton contribution to the n_re equation.
 * Employs the analytical growth rate calculated by RunawayFluid.
 */
using namespace DREAM;
/**
 * Constructor.
 */
ComptonRateTerm::ComptonRateTerm(
    FVM::Grid *g, FVM::UnknownQuantityHandler *uqn,
    RunawayFluid *rf, real_t scaleFactor
) : FVM::DiagonalComplexTerm(g,uqn), REFluid(rf), scaleFactor(scaleFactor) {

    this->AllocateGamma();

    AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD));
    AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_N_TOT));
}

/**
 * Destructor.
 */
ComptonRateTerm::~ComptonRateTerm() {
    DeallocateGamma();
}

/**
 * Allocate memory for the runaway rate.
 */
void ComptonRateTerm::AllocateGamma() {
    this->gamma = new real_t[this->grid->GetNr()];
    this->dGamma = new real_t[this->grid->GetNr()];
}

/**
 * Free memory for the runaway rate.
 */
void ComptonRateTerm::DeallocateGamma() {
    delete [] this->gamma;
    delete [] this->dGamma;
}

/**
 * Method called when the grid has been rebuilt.
 */
bool ComptonRateTerm::GridRebuilt() {
    DeallocateGamma();
    AllocateGamma();
    return FVM::DiagonalComplexTerm::GridRebuilt();
}

void ComptonRateTerm::SetDiffWeights(len_t derivId, len_t nMultiples){
    REFluid->evaluatePartialContributionComptonGrowthRate(dGamma, derivId);
    len_t offset = 0;
    for(len_t n = 0; n<nMultiples; n++){
        for (len_t ir = 0; ir < nr; ir++){
            diffWeights[offset + n1[ir]*(n2[ir]-1) + 0] = scaleFactor*dGamma[ir];
            offset += n1[ir]*n2[ir];
        }
    }
}

void ComptonRateTerm::SetWeights(){
    REFluid->SetComptonRunawayRate(gamma);
    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++){
        weights[offset + n1[ir]*(n2[ir]-1) + 0] = scaleFactor*gamma[ir];
        offset += n1[ir]*n2[ir];
    }
}

