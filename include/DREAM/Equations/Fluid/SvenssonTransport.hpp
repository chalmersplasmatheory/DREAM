#ifndef _DREAM_SVENSSON_TRANSPORT_HPP
#define _DREAM_SVENSSON_TRANSPORT_HPP


#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
//#include "FVM/Interpolator1D.hpp"
//#include "FVM/Interpolator3D.hpp"

namespace DREAM {
    template<typename T>
    class SvenssonTransport : public T {
    protected:
        
        const len_t nr, np, EID;
        const real_t pStar;
        real_t *coeff;
        const real_t *p;
        real_t *integrand;
        
        FVM::UnknownQuantityHandler *unknowns;
        DREAM::RunawayFluid *REFluid;
        FVM::Interpolator3D *interp3d;

        
        void _setcoeff(const len_t, const real_t);
        
        //virtual const real_t *EvaluateIntegrand(len_t)=0;
        virtual void EvaluateIntegrand(len_t)=0;

        real_t GetPBarInv_f(len_t, real_t *dr_pBarInv_f=nullptr);

    public:
        SvenssonTransport<T>(
            FVM::Grid*, real_t,
            FVM::UnknownQuantityHandler*, RunawayFluid*, 
            FVM::Interpolator3D* );
        
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
    class SvenssonTransportDiffusionTerm : public SvenssonTransport<FVM::DiffusionTerm>{
        using SvenssonTransport<FVM::DiffusionTerm>::SvenssonTransport;
        
        // Function for calculating the integrand associated to the
        // diffusion coefficient.
        void EvaluateIntegrand(len_t ir);
    };


    /**
     *   Class for evaluating the integrand associated with the
     *   Gamma-bar (advection like) coefficient in the Svensson
     *   transport. This coefficient only depends on the advection (A)
     *   and diffusion (D) coefficients, so there are two classes
     *   handeling the two input coefficients separately.
     */
    class SvenssonTransportAdvectionTermA : public SvenssonTransport<FVM::AdvectionTerm>{
        using SvenssonTransport<FVM::AdvectionTerm>::SvenssonTransport;

        // Function for calculating the integrand associated to the
        // diffusion coefficient.
        void EvaluateIntegrand(len_t ir);
    };

    /**
     *   Class for evaluating the integrand associated with the
     *   Gamma-bar (advection like) coefficient in the Svensson
     *   transport. This coefficient only depends on the advection (A)
     *   and diffusion (D) coefficients, so there are two classes
     *   handeling the two input coefficients separately.
     */
    class SvenssonTransportAdvectionTermD : public SvenssonTransport<FVM::AdvectionTerm>{
        using SvenssonTransport<FVM::AdvectionTerm>::SvenssonTransport;

        // Function for calculating the integrand associated to the
        // diffusion coefficient.
        void EvaluateIntegrand(len_t ir);
    };
}

#include "SvenssonTransport.tcc"

#endif/*_DREAM_SVENSSON_TRANSPORT_HPP*/


