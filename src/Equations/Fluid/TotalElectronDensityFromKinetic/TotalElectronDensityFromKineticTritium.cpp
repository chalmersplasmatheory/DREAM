#include "DREAM/Equations/Fluid/TotalElectronDensityFromKinetic/TotalElectronDensityFromKineticTritium.hpp"

/**
 * Implementation of an equation term which represents the total
 * number of electrons created by the kinetic Tritium beta decay
 * source
 */ 
using namespace DREAM;

TotalElectronDensityFromKineticTritium::TotalElectronDensityFromKineticTritium(
    FVM::Grid* g, real_t pLower, real_t pUpper, FVM::UnknownQuantityHandler *u, IonHandler *ions, len_t iIon, real_t scaleFactor
) : FVM::EquationTerm(g), pLower(pLower), pUpper(pUpper), scaleFactor(scaleFactor) {

    this->id_nT = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    
    this->indT = ions->GetIndex(iIon, 0);
}

TotalElectronDensityFromKineticTritium::~TotalElectronDensityFromKineticTritium() {
    this->DeallocateWeights();
}

void TotalElectronDensityFromKineticTritium::AllocateWeights() {
    DeallocateWeights(); 

    const len_t N = grid->GetNCells();
    weights = new real_t[N];
    for (len_t i = 0; i < N; i++)
        weights[i] = 0;
}

void TotalElectronDensityFromKineticTritium::DeallocateWeights(){
    if(weights!=nullptr) 
    delete[] weights;
}

void TotalElectronDensityFromKineticTritium::InitializeWeights(){
    AllocateWeights(); 
    SetWeights();
    hasBeenInitialized = false;
}

void TotalElectronDensityFromKineticTritium::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) { 
    if(!hasBeenInitialized){
        InitializeWeights();
        hasBeenInitialized = true;
    } 
}

bool TotalElectronDensityFromKineticTritium::SetJacobianBlock(const len_t , const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(derivId != id_nT)
        return false;
    SetMatrixElements(jac, nullptr);
    return true;
}

void TotalElectronDensityFromKineticTritium::SetMatrixElements(FVM::Matrix *mat, real_t* /*rhs*/) {
    for(len_t iZ0=0; iZ0<2; iZ0++){
        len_t offset = 0;
        for(len_t ir=0; ir<nr; ir++){
            for(len_t i=0; i<n1[ir]; i++){
                for(len_t j=0; j<n2[ir]; j++){
                    len_t ind = offset + n1[ir]*j + i;
                    mat->SetElement(ind, (iZ0 + indT)*nr + ir, weights[ind]);
                }
            }
            offset += n1[ir]*n2[ir];
        }        
    }
}

void TotalElectronDensityFromKineticTritium::SetVectorElements(real_t *vec, const real_t *x) {
    for(len_t iZ0=0; iZ0<2; iZ0++){
        len_t offset = 0; 
        for(len_t ir=0; ir<nr; ir++){
            for(len_t i=0; i<n1[ir]; i++){
                for(len_t j=0; j<n2[ir]; j++){
                    len_t ind = offset + n1[ir]*j + i;
                    vec[ind] += weights[ind]*x[(iZ0 + indT)*nr + ir];
                }
            }
            offset += n1[ir]*n2[ir];
        }
    }
}

void TotalElectronDensityFromKineticTritium::SetWeights() {
    for(len_t i = 0; i<grid->GetNCells(); i++)
        weights[i] = scaleFactor * TritiumSource::EvaluateTotalTritiumNumber(pLower, pUpper);
}