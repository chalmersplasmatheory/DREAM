#ifndef _DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_KINETIC_TERM_HPP
#define _DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_KINETIC_TERM_HPP

#include "FVM/UnknownQuantityHandler.hpp"
//#include "DREAM/Equations/CollisionQuantity.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "FVM/Equation/MomentQuantity.hpp"

namespace DREAM {
    class CollisionalEnergyTransferKineticTerm : public FVM::MomentQuantity {
    private:
        // The collQtyHandler should preferably be the one corresponding to fGrid
        CollisionQuantityHandler *collQtyHandler;
        enum OptionConstants::momentumgrid_type mgtype;
        real_t scaleFactor;
        CollisionQuantity::collqty_settings *collQtySetting;
        virtual void SetDiffIntegrand(len_t derivId) override;

    public:
        CollisionalEnergyTransferKineticTerm(
            FVM::Grid*, FVM::Grid*, len_t, len_t,
            CollisionQuantityHandler*,FVM::UnknownQuantityHandler*,  
            enum OptionConstants::momentumgrid_type, real_t scaleFactor = 1.0,
            real_t pThreshold = 0, pThresholdMode pMode = FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_MC
        );
        virtual ~CollisionalEnergyTransferKineticTerm();

        virtual bool GridRebuilt() override {return false;}
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

    };
}

#endif /*_DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_KINETIC_TERM_HPP*/
