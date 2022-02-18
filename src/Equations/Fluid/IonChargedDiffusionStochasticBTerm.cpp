#include "DREAM/Constants.hpp"
#include "DREAM/Equations/Fluid/IonChargedAdvectionDiffusionTerm.hpp"
#include "DREAM/Equations/Fluid/IonChargedDiffusionStochasticBTerm.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "DREAM/IonHandler.hpp"
//#include "DREAM/NotImplementedException.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

IonChargedDiffusionStochasticBTerm::IonChargedDiffusionStochasticBTerm(FVM::Grid *grid, IonHandler *ions, bool allocCoefficients,
	const len_t iIon, FVM::Interpolator1D* dBOverB, MultiInterpolator1D* DrrHat, FVM::UnknownQuantityHandler *u
	) : IonChargedAdvectionDiffusionTerm<FVM::DiffusionTerm>(grid, 
	ions, iIon, allocCoefficients), dBOverB(dBOverB), DrrHat(DrrHat) {
	
    SetName("IonChargedDiffusionStochasticBTerm");

    this->id_ni = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    AddUnknownForJacobian(u, id_ni);
    
    this->id_Wi = u->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    AddUnknownForJacobian(u, id_Wi);
    
    this->id_Ni = u->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    AddUnknownForJacobian(u, id_Ni);
	
	Allocate();
}

IonChargedDiffusionStochasticBTerm::~IonChargedDiffusionStochasticBTerm(){
	Deallocate();
}

void IonChargedDiffusionStochasticBTerm::Allocate(){
    Deallocate();
    
    len_t nzs=ions->GetNzs();

    const len_t nr = this->grid->GetNr();

	// Derivatives wrt ion (charge state) densities
    this->dDrrdni = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dDrrdni[Z0-1] = new real_t[(nr+1)*nzs];
        
    // Derivatives wrt ion thermal energy (of this species only)
    this->dDrrdWi = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dDrrdWi[Z0-1] = new real_t[(nr+1)];
        
    // Derivatives wrt the total ion density (of this species only)
    this->dDrrdNi = new real_t*[Zion];
    for(len_t Z0=1; Z0<=Zion; Z0++)
        this->dDrrdNi[Z0-1] = new real_t[(nr+1)];
}

void IonChargedDiffusionStochasticBTerm::Deallocate(){
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dDrrdni[Z0-1];
    delete [] this->dDrrdni;
    
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dDrrdWi[Z0-1];
    delete [] this->dDrrdWi;
    
    for(len_t Z0=1; Z0<=Zion; Z0++)
        delete [] this->dDrrdNi[Z0-1];
    delete [] this->dDrrdNi;
}

void IonChargedDiffusionStochasticBTerm::SetCoeffsAllCS(const real_t){
}

void IonChargedDiffusionStochasticBTerm::SetDiffCoeffsAllCS(const real_t){
}

void IonChargedDiffusionStochasticBTerm::SetCoeffs(const len_t Z0){
	if(Z0<1)
		return;
		
	const len_t nr = this->grid->GetNr();
	for(len_t ir=0; ir<nr+1; ir++)
		Drr(ir,0,0)=CoeffsAllCS[Z0-1][ir];
}


void IonChargedDiffusionStochasticBTerm::SetPartialDiffusionTerm(len_t derivId, len_t nMultiples){
	if(derivId==id_ni)
		for(len_t n=0; n<nMultiples; n++)
			for(len_t ir=0; ir<nr+1; ir++)
				dDrr(ir,0,0,n)=dDrrdni[Z0ForPartials-1][ir+nr*n];
				
	if(derivId==id_Wi)
		for(len_t n=0; n<nMultiples; n++)
			if(n==iIon)
				for(len_t ir=0; ir<nr+1; ir++)
					dDrr(ir,0,0,n)=dDrrdWi[Z0ForPartials-1][ir];	
					
	if(derivId==id_Ni)
		for(len_t n=0; n<nMultiples; n++)
			if(n==iIon)
				for(len_t ir=0; ir<nr+1; ir++)
					dDrr(ir,0,0,n)=dDrrdNi[Z0ForPartials-1][ir];			
}




