#ifndef _DREAMTESTS_DREAM_RUNAWAY_FLUID_HPP
#define _DREAMTESTS_DREAM_RUNAWAY_FLUID_HPP

#include <string>
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/ConnorHastie.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
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
        DREAM::FVM::UnknownQuantityHandler *GetUnknownHandlerSingleImpuritySpecies(DREAM::FVM::Grid*, 
            const real_t, const len_t, const len_t, const real_t HYDROGEN_DENSITY=1e20, const real_t T_cold=10);
        DREAM::RunawayFluid *GetRunawayFluid(
            DREAM::CollisionQuantity::collqty_settings*,
            DREAM::CollisionQuantity::collqty_settings*,
            const len_t, const len_t*, const real_t, const real_t, const real_t,
            const len_t nr, enum DREAM::OptionConstants::eqterm_dreicer_mode dm=DREAM::OptionConstants::EQTERM_DREICER_MODE_NONE,
            enum DREAM::OptionConstants::collqty_Eceff_mode em=DREAM::OptionConstants::COLLQTY_ECEFF_MODE_FULL
        );
        DREAM::RunawayFluid *GetRunawayFluidSingleImpuritySpecies(
            DREAM::CollisionQuantity::collqty_settings*, 
            DREAM::CollisionQuantity::collqty_settings*, 
            const real_t IMPURITY_DENSITY, const len_t IMPURITY_Z0, const len_t IMPURITY_Z,
            const real_t B0, 
            enum DREAM::OptionConstants::eqterm_dreicer_mode dm=DREAM::OptionConstants::EQTERM_DREICER_MODE_NONE,
            enum DREAM::OptionConstants::collqty_Eceff_mode em=DREAM::OptionConstants::COLLQTY_ECEFF_MODE_FULL,
            const real_t HYDROGEN_DENSITY=1e20, const real_t T_cold=10
        );
        bool CompareEceffWithTabulated();
        bool CompareGammaAvaWithTabulated();
        bool CompareConnorHastieRateWithTabulated();

        real_t _ConnorHastieFormula(
            const real_t, const real_t, const real_t,
            const real_t, const real_t,
            bool
        );

//        bool CheckConservativity();
        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_DREAM_RUNAWAY_FLUID_HPP*/
