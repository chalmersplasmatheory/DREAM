/**
 * Implementation of the pellet ablation term (takes ablation rate provided by SPIHandler).
 */

#include "DREAM/Equations/Fluid/SPIHeatAbsorbtionTerm.hpp"

using namespace DREAM;

SPIHeatAbsorbtionTerm::SPIHeatAbsorbtionTerm(
    FVM::Grid* g, SPIHandler *SPI, real_t scaleFactor=1.0
): EquationTerm(g),SPI(SPI), scaleFactor(scaleFactor){}

void SPIHeatAbsorbtionTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler* ){
    heatAbsorbtionRate=SPI->GetHeatAbsorbtionRate();
}

bool SPIHeatAbsorbtionTerm::SetJacobianBlock(const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*){
    SPI->setJacobianAdiabaticHeatAbsorbtionRate(jac, derivId, scaleFactor);
    return true;
}


/**
 * Set the linear operator matrix elements corresponding
 * to this term.
 */
void SPIHeatAbsorbtionTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs) {
    this->SetVectorElements(rhs, nullptr);
}

/**
 * Set the non-linear function vector for this term.
 */
void SPIHeatAbsorbtionTerm::SetVectorElements(real_t *vec, const real_t*){
    for(len_t ir=0;ir<nr;ir++){
        vec[ir]+=scaleFactor*heatAbsorbtionRate[ir];
    }
}
