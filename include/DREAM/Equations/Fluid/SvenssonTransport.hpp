#ifndef _DREAM_SVENSSON_TRANSPORT_HPP
#define _DREAM_SVENSSON_TRANSPORT_HPP

#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
//#include "FVM/Interpolator1D.hpp"
//#include "FVM/Interpolator3D.hpp"

namespace DREAM {
    template<typename T>
    class SvenssonTransport : public T {
    private:
        
        const len_t nr, np;
        const real_t pStar;
        const real_t **coeffA, **coeffD, *r, *p1, *p2;
        FVM::UnknownQuantityHandler *unknowns;
        RunawayFluid *REFluid;

        void _setcoeff(const len_t, const len_t, const real_t);
        
        virtual const reat_t *EvaluateIntegrand(len_t)=0;

    public:
        SvenssonTransport<T>(
            FVM::Grid*,
            const len_t, const len_t, const real_t,
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
    void SvenssonTransport<FVM::AdvectionTerm>::_setcoeff(const len_t, const len_t, const real_t);
    template<>
    void SvenssonTransport<FVM::DiffusionTerm>::_setcoeff(const len_t, const len_t, const real_t);

    
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
    class SvenssonTransportDiffusive : public SvenssonTransport<FVM::DiffusionTerm>{
        // public: // ???

        // Add Constructor?
        
        // Function for calculating the integrand associated to the
        // diffusion coefficient.
        const reat_t *EvaluateIntegrand(len_t ir);
    };


    /**
     *   Class for evaluating the integrand associated with the
     *   Gamma-bar (advection like) coefficient in the Svensson
     *   transport. This coefficient only depends on the advection (A)
     *   and diffusion (D) coefficients, so there are two classes
     *   handeling the two input coefficients separately.
     */
    // YYY Work out a more descriptive name
    class SvenssonTransportAdvectiveA : public SvenssonTransport<FVM::AdvectionTerm>{
        // Add constructor??

        // Function for calculating the integrand associated to the
        // diffusion coefficient.
        const reat_t *EvaluateIntegrand(len_t ir);
    };

    /**
     *   Class for evaluating the integrand associated with the
     *   Gamma-bar (advection like) coefficient in the Svensson
     *   transport. This coefficient only depends on the advection (A)
     *   and diffusion (D) coefficients, so there are two classes
     *   handeling the two input coefficients separately.
     */
    // YYY Work out a more descriptive name
    class SvenssonTransportAdvectiveD : public SvenssonTransport<FVM::AdvectionTerm>{


        // Function for calculating the integrand associated to the
        // diffusion coefficient.
        const reat_t *EvaluateIntegrand(len_t ir);
    };
}

#include "SvensonTransport.tcc"

#endif/*_DREAM_SVENSSON_TRANSPORT_HPP*/


