#include "DREAM/Equations/Fluid/TotalElectronDensityFromKinetic/TotalElectronDensityFromKineticAvalancheCH.hpp"

#include "DREAM/Equations/Kinetic/AvalancheSourceCH.hpp"

/**
 * Implementation of an equation term which represents the total
 * number of electrons created by the kinetic Chiu-Harvey source
 */ 
using namespace DREAM;

TotalElectronDensityFromKineticAvalancheCH::TotalElectronDensityFromKineticAvalancheCH(
    FVM::Grid* grid, FVM::Grid* og, real_t pCutoff, FVM::UnknownQuantityHandler *u, 
    AvalancheSourceCH *avaCH, real_t scaleFactor
) : FVM::DiagonalComplexTerm(grid,u,og), pCutoff(pCutoff), avaCH(avaCH), scaleFactor(scaleFactor) {
    id_fhot = this->unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
}
        
void TotalElectronDensityFromKineticAvalancheCH::SetWeights() {
    for(len_t ir=0; ir<grid->GetNr(); ir++){
        weights[ir] = 0;
        len_t np1 = this->operandGrid->GetMomentumGrid(ir)->GetNp1();
        len_t np2 = this->operandGrid->GetMomentumGrid(ir)->GetNp2();
        for(len_t i=0; i<np1; i++){
            for(len_t j=0; j<np2; j++){
                weights[ir] += scaleFactor * this->avaCH->GetSourceFunction(ir, i, j)
                                    * this->operandGrid->GetMomentumGrid(ir)->GetDp1(i) * this->operandGrid->GetMomentumGrid(ir)->GetDp2(j);
            }
        }
    }   
}

bool TotalElectronDensityFromKineticAvalancheCH::AddWeightsJacobian(
    const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t* x
) {
    if (derivId==id_fhot){ 
        for(len_t ir=0; ir<grid->GetNr(); ir++){
            len_t np1 = this->operandGrid->GetMomentumGrid(ir)->GetNp1();
            len_t np2 = this->operandGrid->GetMomentumGrid(ir)->GetNp2();
            for(len_t i=0; i<np1; i++){
                real_t cutFactor = 1.;
                if (this->operandGrid->GetMomentumGrid(ir)->GetP1_f(i+1) < this->pCutoff)
                    cutFactor = 0.;
                else if (this->operandGrid->GetMomentumGrid(ir)->GetP1_f(i) < this->pCutoff)
                    cutFactor = (this->operandGrid->GetMomentumGrid(ir)->GetP1_f(i+1) - pCutoff) / (this->operandGrid->GetMomentumGrid(ir)->GetP1_f(i+1) - this->operandGrid->GetMomentumGrid(ir)->GetP1_f(i));
                for(len_t j=0; j<np2; j++) {
                    real_t derivTerm = scaleFactor * cutFactor * this->avaCH->GetSourceFunctionJacobian(ir, i, j, derivId)
                                            * this->operandGrid->GetMomentumGrid(ir)->GetDp1(i) * this->operandGrid->GetMomentumGrid(ir)->GetDp2(j);
                    jac->SetElement(ir, ir*np2*np1 + j*np1 + i, derivTerm * x[ir]);
                }
            }
        }

        return true;
    }

    return false;
}

/**
* Set the linear operator matrix elements corresponding to this term.
*/
void TotalElectronDensityFromKineticAvalancheCH::SetMatrixElements(FVM::Matrix *mat, real_t*) {
    len_t N = this->grid->GetNCells();
    for (len_t i = 0; i < N; i++)
        mat->SetElement(i, i, weights[i]);
}

/**
* Set function vector for this term.
*/
void TotalElectronDensityFromKineticAvalancheCH::SetVectorElements(real_t *vec, const real_t *x) {
    for (len_t ir = 0; ir < grid->GetNCells(); ir++)
        vec[ir] += weights[ir] * x[ir];
}


bool TotalElectronDensityFromKineticAvalancheCH::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t* x
) {
    bool contributes = false;
    if (derivId == uqtyId) {
        this->SetMatrixElements(jac, nullptr);
        contributes = true;
    }

    contributes |= AddWeightsJacobian(uqtyId, derivId, jac, x);

    return contributes;
}