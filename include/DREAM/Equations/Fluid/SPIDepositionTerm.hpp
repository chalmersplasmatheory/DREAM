#ifndef _DREAM_EQUATION_FLUID_PELLET_HEAT_ABSORBTION_TERM_HPP
#define _DREAM_EQUATION_FLUID_PELLET_HEAT_ABSORBTION_TERM_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class SPIHeatAbsorbtionTerm : public FVM::EquationTerm {

    private:
        SPIHandler *SPI;
        real_t scaleFactor=1.0;

        real_t *heatAbsorbtionRate;

    public:
        SPIHeatAbsorbtionTerm(
            FVM::Grid*, FVM::UnknownQuantityHandler*, SPIHandler*, real_t scaleFactor=1.0
        );
        ~SPIHeatAbsorbtionTerm(){}
        
        //virtual bool GridRebuilt() override;
        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return 1; }   /* XXX TODO XXX */
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_PELLET_HEAT_ABSORBTION_TERM_HPP*/
