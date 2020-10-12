#ifndef _DREAMTESTS_DREAM_MEAN_EXCITATION_ENERGY_HPP
#define _DREAMTESTS_DREAM_MEAN_EXCITATION_ENERGY_HPP

#include <string>
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "UnitTest.hpp"
#include "DREAM/Equations/CollisionQuantity.hpp"

namespace DREAMTESTS::_DREAM {
    class MeanExcitationEnergy : public UnitTest {
    private:
        len_t id_ions;

    public:
        MeanExcitationEnergy(const std::string& s) : UnitTest(s) {}
        virtual ~MeanExcitationEnergy() {}

        DREAM::IonHandler *GetIonHandler(DREAM::FVM::Grid*, DREAM::FVM::UnknownQuantityHandler*, const len_t, const len_t*);
        DREAM::FVM::UnknownQuantityHandler *GetUnknownHandler(DREAM::FVM::Grid*,
            const len_t, const len_t*, const real_t, const real_t);
        real_t *GetMeanExcitationEnergies(DREAM::CollisionQuantity::collqty_settings *cq, const len_t N_IONS, const len_t *Z_IONS,
            const len_t N_SPECIES_TO_TEST, const len_t *Z_TO_TEST, const len_t *Z0_TO_TEST,// real_t *meanExcitationEnergies, 
            const real_t ION_DENSITY_REF, const real_t T_cold, const real_t B0, const len_t nr);
        bool CompareMeanExcitationEnergyWithTabulated();

        virtual bool Run(bool) override;
    };
}

#endif/*_DREAMTESTS_DREAM_MEAN_EXCITATION_ENERGY_HPP*/
