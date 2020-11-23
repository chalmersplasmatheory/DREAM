#ifndef _DREAM_EQUATION_ION_SPI_DEPOSITION_TERM_HPP
#define _DREAM_EQUATION_ION_SPI_DEPOSITION_TERM_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/SPIHandler.hpp"

namespace DREAM {
    class IonSPIDepositionTerm : public IonEquationTerm<FVM::EquationTerm> {
    protected:
        SPIHandler *SPI;
        real_t SPIMolarFraction;
        real_t scaleFactor;
        real_t *weights;
        //real_t *depositionRate;

    public:
        IonSPIDepositionTerm(
            FVM::Grid *g, IonHandler *ihdl, const len_t iIon,
            SPIHandler *SPI, real_t SPIMolarFraction, real_t scaleFactor
        );
        virtual ~IonSPIDepositionTerm();

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 3; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override {return 3;} //One each for ncold, Tcold and rp. NOTE: could be more due to non-locality!

        //virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual void SetCSJacobianBlock(
            const len_t, const len_t, FVM::Matrix*, const real_t*,
            const len_t, const len_t Z0, const len_t rOffset
        ) override;

        virtual void SetCSMatrixElements(
            FVM::Matrix*, real_t*, const len_t, const len_t Z0, const len_t rOffset
        ) override;

        virtual void SetCSVectorElements(
            real_t*, const real_t*, const len_t, const len_t Z0, const len_t rOffset
        ) override;
    };
}

#endif/*_DREAM_EQUATION_ION_SPI_DEPOSITION_TERM_HPP*/

