#ifndef _DREAMTESTS_DREAM_RUNAWAY_FLUID_HPP
#define _DREAMTESTS_DREAM_RUNAWAY_FLUID_HPP

#include <string>
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "UnitTest.hpp"

namespace DREAMTESTS::_DREAM {
    class RunawayFluid : public UnitTest {
    private:
        len_t id_ions;

    public:
        RunawayFluid(const std::string& s) : UnitTest(s) {}

        DREAM::IonHandler *GetIonHandler(DREAM::FVM::Grid*, DREAM::FVM::UnknownQuantityHandler*, const len_t, const len_t*);
        DREAM::FVM::UnknownQuantityHandler *GetUnknownHandler(DREAM::FVM::Grid*,
            const len_t, const len_t*, const real_t, const real_t);
        DREAM::RunawayFluid *GetRunawayFluid(DREAM::CollisionQuantity::collqty_settings *cq, const len_t, const len_t*,const real_t, const real_t, const real_t, const len_t nr);
        bool CompareEceffWithTabulated();
        bool CompareGammaAvaWithTabulated();

//        bool CheckConservativity();
        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_DREAM_RUNAWAY_FLUID_HPP*/
