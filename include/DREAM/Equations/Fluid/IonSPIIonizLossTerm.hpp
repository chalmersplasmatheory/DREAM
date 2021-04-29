#ifndef _DREAM_EQUATION_ION_SPI_IONIZ_LOSS_TERM_HPP
#define _DREAM_EQUATION_ION_SPI_IONIZ_LOSS_TERM_HPP

/**
 * Implementation of the ionization losses corresponding to the SPI source, based on the data provided by the SPIHandler.
 * This term makes it possible to maintain energy conservation while still being allowed to add the ablated material to a charged state,
 * for example since the material is not really confined until it is ionized, and it can also be helpfull to avoid having to resolv
 * the very rapid ionization through the first charge states
 *
 * Note that this equation is applied to a single _ion species_,
 * (and to all its charge states).
 */

#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/Equations/Fluid/IonSPIDepositionTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/SPIHandler.hpp"
#include "DREAM/Constants.hpp"
#include "DREAM/NIST.hpp"

namespace DREAM {
    class IonSPIIonizLossTerm : public IonSPIDepositionTerm {
    private:
    real_t *EIonizTot;
    public:
        IonSPIIonizLossTerm(
            FVM::Grid *g, IonHandler *ihdl, const len_t iIon,ADAS *adas, FVM::UnknownQuantityHandler *unknowns,bool addFluidIonization, bool addFluidJacobian,
            SPIHandler *SPI, const real_t *SPIMolarFraction, len_t offset, real_t scaleFactor, NIST *nist, bool isAbl = false,
            OptionConstants::eqterm_spi_abl_ioniz_mode spi_abl_ioniz_mode = OptionConstants::EQTERM_SPI_ABL_IONIZ_MODE_SELF_CONSISTENT
        ) : IonSPIDepositionTerm(g, ihdl, iIon, adas, unknowns, addFluidIonization, addFluidJacobian, SPI, SPIMolarFraction, offset, scaleFactor, isAbl, spi_abl_ioniz_mode) {
        	EIonizTot = new real_t[Zion+1];
            for (len_t i=0;i<Zion+1;i++)
                EIonizTot[i]=(nist->GetBindingEnergy(Zion,0)-nist->GetBindingEnergy(Zion,i))*Constants::ec;
        }
        
        virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler* unknowns) override{
        	IonSPIDepositionTerm::Rebuild(t,dt,unknowns);
        	len_t Nr=grid->GetNr();
        	for(len_t ir=0;ir<Nr;ir++)
        		for(len_t iZ=0;iZ<Zion+1;iZ++)
        			weights[ir*(Zion+1)+iZ]*=EIonizTot[iZ];
		}


        virtual bool SetCSJacobianBlock(
            const len_t uqtyId, const len_t derivId, FVM::Matrix* jac, const real_t* x,
            const len_t iIon, const len_t Z0, const len_t
        ) override {
            IonSPIDepositionTerm::SetCSJacobianBlock(uqtyId,derivId,jac,x,iIon,Z0,0);
        }

        virtual void SetCSMatrixElements(
            FVM::Matrix* mat, real_t* rhs, const len_t iIon, const len_t Z0, const len_t
        ) override {
            IonSPIDepositionTerm::SetCSMatrixElements(mat,rhs,iIon,Z0,0);
        }

        virtual void SetCSVectorElements(
            real_t* vec, const real_t* x, const len_t iIon, const len_t Z0, const len_t
        ) override {
            IonSPIDepositionTerm::SetCSVectorElements(vec,x,iIon,Z0,0);
        }
    };
}

#endif/*_DREAM_EQUATION_ION_SPI_IONIZ_LOSS_TERM_HPP*/

