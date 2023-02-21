#ifndef _DREAM_EQUATIONS_BOOTSTRAP_CURRENT_HPP
#define _DREAM_EQUATIONS_BOOTSTRAP_CURRENT_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CoulombLogarithm.hpp"
#include "DREAM/Equations/CollisionQuantity.hpp"
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"

namespace DREAM {
    class BootstrapEquationTerm {
    private:
        FVM::RadialGrid *rGrid;
        FVM::UnknownQuantityHandler *unknowns;
        IonHandler *ions;

        CoulombLogarithm lnLambda;
        struct CollisionQuantity::collqty_settings *lnLEE_settings;
        struct CollisionQuantity::collqty_settings *lnLII_settings;

        real_t epsilon; // used for numerical differentiation

        len_t
            id_jtot,
            id_ncold,
            id_ni,
            id_Tcold,
            id_Wi;

        real_t
            *jtot,
            *niMain,
            *WiMain;

        real_t *ft=nullptr;              // fraction of trapped particles
        real_t *qR0=nullptr;             // safety factor normalised to R0

        real_t
            evaluateElectronCollisionFrequency(len_t, real_t, real_t, real_t),
            evaluateElectronCollisionFrequency(len_t),
            evaluateIonCollisionFrequency(len_t, real_t, real_t, real_t),
            evaluateIonCollisionFrequency(len_t);


    protected:
        bool includeIonTemperatures;

        real_t
            *ncold,
            *Tcold,
            *ni,
            *Wi,
            *p,
            *n;


        real_t
            *constantPrefactor=nullptr,
            *coefficientL31=nullptr;

        /**
         * Coefficients L31, L32 and alpha (as defined in Redl et al. 2021)
         */
        real_t
            evaluateCoefficientL31(len_t, real_t, real_t, real_t, real_t),
            evaluateCoefficientL31(len_t),
            evaluateCoefficientL32(len_t, real_t, real_t, real_t, real_t),
            evaluateCoefficientL32(len_t),
            evaluateCoefficientAlpha(len_t, real_t, real_t),
            evaluateCoefficientAlpha(len_t);

        /**
         * Partial derivatives of L31, L32 and alpha.
         */
        real_t
            evaluateNumericalDerivative(len_t, len_t, len_t, real_t, std::function<real_t(len_t, real_t, real_t, real_t)>),
            evaluatePartialCoefficientL31(len_t, len_t, len_t),
            evaluatePartialCoefficientL32(len_t, len_t, len_t),
            evaluatePartialCoefficientAlpha(len_t, len_t, len_t);


    public:
        BootstrapCurrent(FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*);
        ~BootstrapCurrent();

        void
            AllocateQuantities(),
            DeallocateQuantities(),
            Rebuild();
    };
}



#endif /*_DREAM_EQUATIONS_BOOTSTRAP_CURRENT_HPP*/
