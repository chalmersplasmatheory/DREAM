#ifndef _DREAM_EQUATION_FLUID_PELLET_HEAT_ABSORBTION_TERM_HPP
#define _DREAM_EQUATION_FLUID_PELLET_HEAT_ABSORBTION_TERM_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/Equations/SPIHandler.hpp"

namespace DREAM {
    class SPIHeatAbsorbtionTerm : public FVM::EquationTerm {

    private:
        SPIHandler *SPI;
        real_t scaleFactor;

        real_t *heatAbsorbtionRate;

    public:
        SPIHeatAbsorbtionTerm(
            FVM::Grid*, SPIHandler*, real_t scaleFactor
        );
        ~SPIHeatAbsorbtionTerm(){}
        
        //virtual bool GridRebuilt() override;
        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return 3; }   /*1 for ncold, 1 for Tcold and 1 for rp. NOTE: could be more due to non-locality! */
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_PELLET_HEAT_ABSORBTION_TERM_HPP*/
