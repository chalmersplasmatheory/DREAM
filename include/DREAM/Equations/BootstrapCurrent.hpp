#ifndef _DREAM_EQUATIONS_BOOTSTRAP_CURRENT_HPP
#define _DREAM_EQUATIONS_BOOTSTRAP_CURRENT_HPP

namespace DREAM { class BootstrapCurrent; }

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"
#include "DREAM/Equations/CollisionQuantity.hpp"
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"

namespace DREAM {
    class BootstrapCurrent {
    private:
        FVM::RadialGrid *rGrid;
        FVM::UnknownQuantityHandler *unknowns;
        IonHandler *ions;

        CoulombLogarithm *lnLambda;
        struct CollisionQuantity::collqty_settings *lnLEE_settings;
        struct CollisionQuantity::collqty_settings *lnLII_settings;

        // used for numerical differentiation
        real_t epsilon_central;     // central difference (nu)
        real_t epsilon_forward;     // forward diffeerence (Zeff)

        len_t iZMain;   // main ion species
        len_t nr;

        len_t id_jtot;
        len_t id_ncold;
        len_t id_ions;
        len_t id_Ni;
        len_t id_Tcold;
        len_t id_Wi;

        real_t *jtot;
        real_t *NiMain=nullptr;
        real_t *WiMain=nullptr;

        real_t *ft=nullptr;             // fraction of trapped particles
        real_t *qR0=nullptr;            // safety factor multiplied by R0
        real_t *eps=nullptr;            // safety factor multiplied by R0

        real_t *constantPrefactor=nullptr;

        real_t *coefficientL31=nullptr;
        real_t *coefficientL32=nullptr;
        real_t *coefficientAlpha=nullptr;

        // bool qFromCurrent = false;

        static real_t evaluateCoefficientL31(real_t, real_t, real_t);
        static real_t evaluateCoefficientL32(real_t, real_t, real_t);
        static real_t evaluateCoefficientAlpha(real_t, real_t, real_t);

    public:
        bool includeIonTemperatures;

        real_t *ncold;
        real_t *Tcold;
        real_t *Ni;
        real_t *Wi;
        real_t *p;
        real_t *n;

        BootstrapCurrent(FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*, CoulombLogarithm*, enum OptionConstants::eqterm_bootstrap_mode bootstrap_mode);
        ~BootstrapCurrent();

        void AllocateQuantities();
        void DeallocateQuantities();
        void Rebuild();

        real_t evaluateElectronCollisionFrequency(len_t);
        real_t evaluateIonCollisionFrequency(len_t);

        const real_t getConstantPrefactor(len_t ir) const
            { return constantPrefactor[ir]; }
        const real_t* getConstantPrefactor() const
            { return constantPrefactor; }

        const real_t getCoefficientL31(len_t ir) const
            { return coefficientL31[ir]; }
        const real_t* getCoefficientL31() const
            { return coefficientL31; }

        const real_t getCoefficientL32(len_t ir) const
            { return coefficientL32[ir]; }
        const real_t* getCoefficientL32() const
            { return coefficientL32; }

        const real_t getCoefficientAlpha(len_t ir) const
            { return coefficientAlpha[ir]; }
        const real_t* getCoefficientAlpha() const
            { return coefficientAlpha; }

        /**
         * Partial derivatives of L31, L32 and alpha.
         */
        real_t evaluateNumericalDerivative(len_t, len_t, len_t, std::function<real_t(real_t, real_t, real_t)>);
        real_t evaluatePartialCoefficientL31(len_t, len_t, len_t);
        real_t evaluatePartialCoefficientL32(len_t, len_t, len_t);
        real_t evaluatePartialCoefficientAlpha(len_t, len_t, len_t, len_t);

    };
}



#endif /*_DREAM_EQUATIONS_BOOTSTRAP_CURRENT_HPP*/
