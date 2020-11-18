#ifndef _DREAM_SVENSSON_TRANSPORT_HPP
#define _DREAM_SVENSSON_TRANSPORT_HPP


#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/Interpolator3D.hpp"

namespace DREAM {
    template<typename T>
    class SvenssonTransport : public T {
    protected:
        
        const len_t nr_f, nt, nr, np1, np2, np, nxi, EID;
        const real_t pStar;
        const real_t *t, *r, *p1, *p2, *xi;
        real_t *p;
        
        real_t *coeffTRXiP,     // Size nt*nr_f*nxi*np
            *coeffRP;           // Size nr_f*np
        real_t **coeff4dInput;  // Size nt-by-(nr*np2*np1)
        real_t *integrand;      // Size np

        // Type of momentum grid used for the input data
        enum FVM::Interpolator3D::momentumgrid_type inputMomentumGridType;
        enum FVM::Interpolator3D::interp_method inputInterp3dMethod;
        
        FVM::UnknownQuantityHandler *unknowns;
        DREAM::RunawayFluid *REFluid;
        enum FVM::Interpolator1D::interp_method timeInterpMethod;
        DREAM::FVM::Interpolator1D *interpTCoeff = nullptr;
        
        
        void _setcoeff(const len_t, const real_t);

        // Function for setting (const) np in the init list
        const len_t CountNp(const len_t, const real_t, const real_t*);
        // Function for setting (const) nxi in the init list
        const len_t CountNxi(const len_t np2In){ return np2In; }
            
        void SetMomentumCoordinate();
            
        void InterpolateCoefficient();

        void xiAverage(const real_t*);

        virtual void EvaluateIntegrand(len_t)=0;

        real_t GetPBarInv_f(len_t, real_t *dr_pBarInv_f=nullptr);

        real_t EvalOnFluxGrid(len_t, const real_t*);

    public:
        // Constructor
        SvenssonTransport<T>( 
            FVM::Grid*, real_t,
            FVM::UnknownQuantityHandler*, RunawayFluid*, 
            struct dream_4d_data*,
            enum FVM::Interpolator1D::interp_method timeInterpMethod = FVM::Interpolator1D::interp_method::INTERP_NEAREST );
        
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
     *  Classes for evaluating the integrand associated with the
     *  diffusion- and advection-like terms in the Svensson transport
     *  model. Each class contains a function for calculating the
     *  integrand associated with each term and each coefficient,
     *  respectively.

     *  
     *  The diffusion term only depends on the diffusion (D)
     *  coefficient. Whlie the advection term depends on both the
     *  advection (A) and diffusion (D) coefficients, so there are two
     *  classes handling the two input coefficients separately.
     */
    class SvenssonTransportDiffusionTerm : public SvenssonTransport<FVM::DiffusionTerm>{
        using SvenssonTransport<FVM::DiffusionTerm>::SvenssonTransport;        
        void EvaluateIntegrand(len_t ir);
    };

    class SvenssonTransportAdvectionTermA : public SvenssonTransport<FVM::AdvectionTerm>{
        using SvenssonTransport<FVM::AdvectionTerm>::SvenssonTransport;
        void EvaluateIntegrand(len_t ir);
    };

    class SvenssonTransportAdvectionTermD : public SvenssonTransport<FVM::AdvectionTerm>{
        using SvenssonTransport<FVM::AdvectionTerm>::SvenssonTransport;
        void EvaluateIntegrand(len_t ir);
    };
}

#include "SvenssonTransport.tcc"

#endif/*_DREAM_SVENSSON_TRANSPORT_HPP*/


