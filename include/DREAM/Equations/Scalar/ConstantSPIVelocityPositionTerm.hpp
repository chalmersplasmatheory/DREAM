#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"

#ifndef _DREAM_FVM_EQUATION_CONSTANT_SPI_VELOCITY_POSITION_TERM_HPP
#define _DREAM_FVM_EQUATION_CONSTANT_SPI_VELOCITY_POSITION_TERM_HPP

namespace DREAM {
    class ConstantSPIVelocityPositionTerm : public FVM::PredeterminedParameter {
    private:
        const real_t *initial_data=nullptr;
        real_t *currentData=nullptr;

        len_t id_vp;
        len_t nShard;
        real_t *vp=nullptr;
        
        const real_t *t_delay = nullptr;

    public:
        ConstantSPIVelocityPositionTerm(FVM::Grid *g, FVM::UnknownQuantityHandler* u, const real_t *xp,  const real_t *t_delay = nullptr);
        ~ConstantSPIVelocityPositionTerm();


        const real_t *GetData() { return currentData; }
        void SetData(const real_t* xp_init, bool copy=true);
        void SetTimeDelay(const real_t *t_delay, bool copy=true);
        virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler* u) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override{}
        virtual void SetVectorElements(real_t*, const real_t*) override{}
    };
}

#endif/*_DREAM_FVM_EQUATION_CONSTANT_SPI_VELOCITY_TERM_HPP*/
