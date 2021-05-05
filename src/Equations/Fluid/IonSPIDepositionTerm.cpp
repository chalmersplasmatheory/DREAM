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
    FVM::Grid *g, IonHandler *ihdl, const len_t iIon, ADAS *adas, FVM::UnknownQuantityHandler *unknowns,bool addFluidIonization, bool addFluidJacobian,
    SPIHandler *SPI, const real_t *SPIMolarFraction, len_t offset, real_t scaleFactor = 1.0, bool isAbl = false,
    OptionConstants::eqterm_spi_abl_ioniz_mode spi_abl_ioniz_mode = OptionConstants::EQTERM_SPI_ABL_IONIZ_MODE_SINGLY_IONIZED
) : IonRateEquation(g, ihdl, iIon, adas, unknowns, addFluidIonization, addFluidJacobian, isAbl ), SPI(SPI), scaleFactor(scaleFactor){

	len_t Nr = this->grid->GetNr();
    weights = new real_t[(Zion+1)*Nr];
    weightsCS = new real_t[Nr];
    //for (len_t i=0;i<(Zion+1)*Nr;i++)
    //    weights[i]=1;
        
    len_t nShard = SPI->GetNShard();
    this->SPIMolarFraction = new real_t[nShard];
    for(len_t ip=0;ip<nShard;ip++)
    	this->SPIMolarFraction[ip] = SPIMolarFraction[offset+ip];
    	
    this->isAbl=isAbl;
    this->spi_abl_ioniz_mode=spi_abl_ioniz_mode;
}

/**
 * Destructor.
 */
IonSPIDepositionTerm::~IonSPIDepositionTerm() {
	delete [] weights;
	delete [] SPIMolarFraction;
}

/**
 * Rebuild rate coefficients.
 */
void IonSPIDepositionTerm::Rebuild(
    const real_t t, const real_t dt, FVM::UnknownQuantityHandler* unknowns
) {

	len_t Nr=this->grid->GetNr();
	if(this->isAbl && this->spi_abl_ioniz_mode==OptionConstants::EQTERM_SPI_ABL_IONIZ_MODE_SINGLY_IONIZED){
		for(len_t ir=0;ir<Nr;ir++)
			for(len_t iZ=0;iZ<Zion+1;iZ++)
				if(iZ==1)
					weights[ir*(Zion+1)+iZ]=1;
				else
					weights[ir*(Zion+1)+iZ]=0;
	
	}else{
		IonRateEquation::Rebuild(t,dt,unknowns);
		real_t sum;
		for(len_t ir=0;ir<Nr;ir++){
			weights[ir*(Zion+1)]=1;
			sum=1;
			for(len_t iZ=0;iZ<Zion;iZ++){
				weights[ir*(Zion+1)+iZ+1]=weights[ir*(Zion+1)+iZ]*Ion[iZ][ir]/Rec[iZ+1][ir];
				sum+=weights[ir*(Zion+1)+iZ+1];
			}
			for(len_t iZ=0;iZ<Zion+1;iZ++)
				weights[ir*(Zion+1)+iZ]*=1/sum;
		}
	}
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
bool IonSPIDepositionTerm::SetCSJacobianBlock(
    const len_t, const len_t derivId, FVM::Matrix *jac,
    const real_t* ,
    const len_t, const len_t Z0, const len_t rOffset
) {

	len_t Nr=this->grid->GetNr();
	for(len_t ir=0;ir<Nr;ir++)
		weightsCS[ir]=scaleFactor*weights[ir*(Zion+1)+Z0];
		
	SPI->evaluatePartialContributionDepositionRate(jac,derivId, weightsCS, SPIMolarFraction, rOffset);
    //if(Z0==Zion){//All deposited material of this ion species is added to the first charge state
    //    SPI->evaluatePartialContributionDepositionRate(jac,derivId, weights[Z0]*scaleFactor, SPIMolarFraction, rOffset);
    //}
    
    return true;
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

    //if(Z0==Zion){
        real_t *depositionRate = SPI->CalculateDepositionRate(SPIMolarFraction);
        const len_t nr = this->grid->GetNr();
        for(len_t ir=0;ir<nr;ir++){
            vec[rOffset+ir]+=scaleFactor*weights[ir*(Zion+1)+Z0]*depositionRate[ir];
        }
    //}
}
