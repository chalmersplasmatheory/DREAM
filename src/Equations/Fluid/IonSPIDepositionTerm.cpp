/**
 * Implementation of the SPI source, based on the data provided by the SPIHandler.
 * The deposited material is added to the singly ionized charge state
 *
 * Note that this equation is applied to a single _ion species_,
 * (and to all its charge states).
 */
#include <iostream>

#include "DREAM/Equations/Fluid/IonSPIDepositionTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
IonSPIDepositionTerm::IonSPIDepositionTerm(
    FVM::Grid *g, IonHandler *ihdl, const len_t iIon,
    SPIHandler *SPI, real_t SPIMolarFraction, real_t scaleFactor = 1.0
) : IonEquationTerm<FVM::EquationTerm>(g, ihdl, iIon), SPI(SPI), SPIMolarFraction(SPIMolarFraction),scaleFactor(scaleFactor){}

/**
 * Destructor.
 */
IonSPIDepositionTerm::~IonSPIDepositionTerm() {}




/**
 * Rebuild rate coefficients.
 */
void IonSPIDepositionTerm::Rebuild(
    const real_t, const real_t, FVM::UnknownQuantityHandler*
) {

    /* //cout<<depositionRate<<endl;
    cout<<"Updating deposition rate"<<endl;
    real_t *depositionRate=SPI->GetDepositionRate();
    //depositionRate=SPI->GetDepositionRate();
    cout<<"Deposition rate updated"<<endl;*/
}


/**
 * Build block of Jacobian matrix for the given charge state.
 *
 * derivId: ID of unknown quantity with respect to which differentiation
 *          should be carried out.
 * uqtyId:  ID of unknown quantity to differentiate.
 * jac:     Jacobian matrix to build.
 * x:       Current value of the unknown quantity.
 * iIon:    Index of ion to build jacobian for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonSPIDepositionTerm::SetCSJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac,
    const real_t* ,
    const len_t, const len_t Z0, const len_t rOffset
) {
    if(Z0==1){//All deposited material of this ion species is added to the first charge state
        SPI->evaluatePartialContributionDepositionRate(jac,derivId, scaleFactor, SPIMolarFraction, rOffset);
    }
}

/**
 * Build linear operator matrix for this equation.
 *
 * mat:     Linear operator matrix.
 * rhs:     Vector representing the equation right-hand-side.
 * iIon:    Index of ion to build matrix for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonSPIDepositionTerm::SetCSMatrixElements(
    FVM::Matrix*, real_t* rhs, const len_t iIon, const len_t Z0, const len_t rOffset
) {
    this->SetCSVectorElements(rhs, nullptr, iIon, Z0, rOffset);
}

/**
 * Build function vector.
 *
 * vec:     Function vector to set elements of.
 * nions:   Ion densities.
 * iIon:    Index of ion to build matrix for.
 * Z0:      Ion charge state.
 * rOffset: Offset in matrix block to set elements of.
 */
void IonSPIDepositionTerm::SetCSVectorElements(
    real_t *vec, const real_t*,
    const len_t, const len_t Z0, const len_t rOffset
) {

    if(Z0==1){
        real_t *depositionRate = SPI->GetDepositionRate();
        const len_t nr = this->grid->GetNr();
        for(len_t ir=0;ir<nr;ir++){
            vec[rOffset+ir]+=scaleFactor*SPIMolarFraction*depositionRate[ir];
        }
    }
}

