/**
 * Implementation of the pellet ablation term (takes ablation rate provided by SPIHandler).
 */

#include "DREAM/Equations/Fluid/SPIAblationTerm.hpp"

using namespace DREAM;

SPIAblationTerm::SPIAblationTerm(
    FVM::Grid* g, FVM::UnknownQuantityHandler *uqn, SPIHandler *SPI, real_t scaleFactor=1.0
): EquationTerm(g),SPI(SPI), scaleFactor(scaleFactor){

    this->id_Yp=uqn->GetUnknownID(OptionConstants::UQTY_Y_P);

}

void SPIAblationTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler* uqn){
    nShard=uqn->GetUnknown(id_Yp)->NumberOfMultiples();
    Ypdot=SPI->GetYpdot();
}

bool SPIAblationTerm::SetJacobianBlock(const len_t, const len_t derivId, FVM::Matrix *jac, const real_t*){
    SPI->setJacobianYpdot(jac, derivId, scaleFactor);
    
    return true;
}


/**
 * Set the linear operator matrix elements corresponding
 * to this term.
 */
void SPIAblationTerm::SetMatrixElements(FVM::Matrix*, real_t *rhs) {
    this->SetVectorElements(rhs, nullptr);
}

/**
 * Set the non-linear function vector for this term.
 */
void SPIAblationTerm::SetVectorElements(real_t *vec, const real_t*){
    for(len_t ip=0;ip<nShard;ip++){
        vec[ip]+=scaleFactor*Ypdot[ip];
    }
}
