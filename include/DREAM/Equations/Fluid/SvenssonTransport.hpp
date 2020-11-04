#ifndef _DREAM_SVENSSON_TRANSPORT_HPP
#define _DREAM_SVENSSON_TRANSPORT_HPP

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
//#include "FVM/Interpolator1D.hpp"
//#include "FVM/Interpolator3D.hpp"

namespace DREAM {
    template<typename T>
    class SvenssonTransport : public T {
    protected:
        
        const len_t nr, np;
        len_t EID;
        const real_t pStar;
        const real_t **coeffA, **coeffD, *r, *p;
        FVM::UnknownQuantityHandler *unknowns;
        DREAM::RunawayFluid *REFluid;

        real_t *integrand;
        
        void _setcoeff(const len_t, const real_t);
        
        virtual const real_t *EvaluateIntegrand(len_t)=0;

        real_t GetPBarInv_f(len_t, real_t *dr_pBarInv_f=nullptr);

    public:
        SvenssonTransport<T>(
            FVM::Grid*,
            const len_t, const len_t,
            const real_t,
            const real_t**, const real_t**, const real_t*, const real_t*,
            FVM::UnknownQuantityHandler*,
            RunawayFluid*,
            bool allocCoefficients=false
        );
        virtual ~SvenssonTransport<T>();

        const real_t *GetCoefficient(const len_t ir);
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };

    template<>
    void SvenssonTransport<FVM::AdvectionTerm>::_setcoeff(const len_t, const real_t);
    template<>
    void SvenssonTransport<FVM::DiffusionTerm>::_setcoeff(const len_t, const real_t);

    
    // Typedefs
    typedef SvenssonTransport<FVM::AdvectionTerm> SvenssonTransportAdvective;
    typedef SvenssonTransport<FVM::DiffusionTerm> SvenssonTransportDiffusive;



    
    /**
     *   Class for evaluating the integrand associated with the
     *   Gamma-tilde (diffusion like) coefficient in the Svensson
     *   transport. This coefficient only depends on the diffusion
     *   coefficient.
     */
    // YYY Work out a more descriptive name!
    class SvenssonTransportDiffusionTerm : public SvenssonTransport<FVM::DiffusionTerm>{
        // Function for calculating the integrand associated to the
        // diffusion coefficient.
        const real_t *EvaluateIntegrand(len_t ir);
    };


    /**
     *   Class for evaluating the integrand associated with the
     *   Gamma-bar (advection like) coefficient in the Svensson
     *   transport. This coefficient only depends on the advection (A)
     *   and diffusion (D) coefficients, so there are two classes
     *   handeling the two input coefficients separately.
     */
    // YYY Work out a more descriptive name
    class SvenssonTransportAdvectionTermA : public SvenssonTransport<FVM::AdvectionTerm>{
        // Function for calculating the integrand associated to the
        // diffusion coefficient.
        const real_t *EvaluateIntegrand(len_t ir);
    };

    /**
     *   Class for evaluating the integrand associated with the
     *   Gamma-bar (advection like) coefficient in the Svensson
     *   transport. This coefficient only depends on the advection (A)
     *   and diffusion (D) coefficients, so there are two classes
     *   handeling the two input coefficients separately.
     */
    // YYY Work out a more descriptive name
    class SvenssonTransportAdvectionTermD : public SvenssonTransport<FVM::AdvectionTerm>{
        // Function for calculating the integrand associated to the
        // diffusion coefficient.
        const real_t *EvaluateIntegrand(len_t ir);
    };
}

//#include "DREAM/Equations/Fluid/SvensonTransport.tcc"
#include "SvenssonTransport.tcc"

#endif/*_DREAM_SVENSSON_TRANSPORT_HPP*/


