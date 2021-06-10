#ifndef _DREAM_EQUATION_FLUID_PELLET_ABLATION_TERM_HPP
#define _DREAM_EQUATION_FLUID_PELLET_ABLATION_TERM_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/Equations/SPIHandler.hpp"

namespace DREAM {
    class SPIAblationTerm : public FVM::EquationTerm {

    private:
        SPIHandler *SPI;
        real_t scaleFactor;

        real_t *Ypdot;
        len_t id_Yp;
        len_t nShard;

    public:
        SPIAblationTerm(
            FVM::Grid*, FVM::UnknownQuantityHandler*,SPIHandler*, real_t scaleFactor
        );
        ~SPIAblationTerm(){}
        
        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return 3; }   /* One each for Tcold, ncold and rp */
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_PELLET_ABLATION_TERM_HPP*/
