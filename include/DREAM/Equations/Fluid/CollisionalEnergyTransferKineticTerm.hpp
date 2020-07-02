#ifndef _DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_KINETIC_TERM_HPP
#define _DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_KINETIC_TERM_HPP

#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/CollisionQuantity.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "FVM/Equation/MomentQuantity.hpp"

namespace DREAM {
    class CollisionalEnergyTransferKineticTerm : public FVM::MomentQuantity {
    private:
        // The collQtyHandler should preferably be the one corresponding to fGrid
        CollisionQuantityHandler *collQtyHandler;
        CollisionQuantity::collqty_settings *collQtySetting;
    public:
        CollisionalEnergyTransferKineticTerm(FVM::Grid*, FVM::Grid*, len_t, len_t,
            CollisionQuantityHandler*);
        virtual ~CollisionalEnergyTransferKineticTerm();

        virtual bool GridRebuilt() override {return false;}
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

    };
}

#endif /*_DREAM_EQUATION_FLUID_COLLISIONAL_ENERGY_TRANSFER_KINETIC_TERM_HPP*/
