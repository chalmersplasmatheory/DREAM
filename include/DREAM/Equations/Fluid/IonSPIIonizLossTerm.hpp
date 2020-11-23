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
    public:
        IonSPIIonizLossTerm(
            FVM::Grid *g, IonHandler *ihdl, const len_t iIon,
            SPIHandler *SPI, real_t SPIMolarFraction, real_t scaleFactor, NIST *nist
        ) : IonSPIDepositionTerm(g, ihdl, iIon, SPI, SPIMolarFraction, scaleFactor) {
            for (len_t i=0;i<Zion+1;i++)
                weights[i]=(nist->GetBindingEnergy(Zion,0)-nist->GetBindingEnergy(Zion,i))*Constants::ec;
        }

        virtual void SetCSJacobianBlock(
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

