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
        // Diffusion coefficient data, prescribed in time and
        // on the given computational grid.
        FVM::Interpolator1D *prescribedCoeff=nullptr;

        len_t nt, nr, np1, np2;
        const real_t **coeffA, **coeffD, *t, *r, *p1, *p2;
        enum FVM::Interpolator3D::momentumgrid_type
            momtype,    // Type of momentum grid used for 'coeff'
            gridtype;   // Type of momentum grid used for 'grid'
        enum FVM::Interpolator3D::interp_method interpmethod = FVM::Interpolator3D::INTERP_LINEAR;

        real_t **interpolateddata=nullptr;

        void _setcoeff(const len_t, const len_t, const real_t);

        virtual const reat_t* EvaluateIntegrand(len_t)=0;

    public:
        SvenssonTransport<T>(
            FVM::Grid*,
            const len_t, const len_t, const len_t, const len_t,
            const real_t**, const real_t**, const real_t*, const real_t*,
            const real_t*, const real_t*, enum FVM::Interpolator3D::momentumgrid_type,
            enum FVM::Interpolator3D::momentumgrid_type,
            enum FVM::Interpolator3D::interp_method interpmethod=FVM::Interpolator3D::INTERP_LINEAR,
            bool allocCoefficients=false
        );
        virtual ~SvenssonTransport<T>();

        const real_t *GetCoefficient(const len_t ir);

      //void InterpolateCoefficient();

      //virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };

    template<>
    void DREAM::SvenssonTransport<DREAM::FVM::AdvectionTerm>::_setcoeff(const len_t, const len_t, const real_t);
    template<>
    void DREAM::SvenssonTransport<DREAM::FVM::DiffusionTerm>::_setcoeff(const len_t, const len_t, const real_t);

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
    class SvenssonTransportGammaTilde : public SvenssonTransport<FVM::DiffusionTerm>{
	// Function for calculating the integrand associated to the
	// diffusion coefficient
	const reat_t* EvaluateIntegrand(len_t ir){
	    // The E field (array over nr) are taken from:
	    real_t *E = this->unknowns->GetUnknownData(this->EID);
	    // The rest of the quantities needed can be foun from:
	    // this->REFluid->Get...()
	    for( len_t i=0; i < this->np; i++){
		pBar = E[ir];
	    }
	    return ...; // An array (in p) with the value of the integrands
	}
    };

    /**
     *   Class for evaluating the integrand associated with the
     *   Gamma-bar (advection like) coefficient in the Svensson
     *   transport. This coefficient only depends on the advection (A)
     *   and diffusion (D) coefficients, so there are two classes
     *   handeling the two input coefficients separately.
     */
    // YYY Work out a more descriptive name
    class SvenssonTransportGammaBarA : public SvenssonTransport<FVM::AdvectionTerm>{
	// Function for calculating the integrand associated to the
	// advection coefficient
	const reat_t* EvaluateIntegrand(len_t ir){
	    // The E field (array over nr) are taken from:
	    real_t *E = this->unknowns->GetUnknownData(this->EID);
	    // The rest of the quantities needed can be foun from:
	    // this->REFluid->Get...()
	    for( len_t i=0; i < this->np; i++){
		pBar = E[ir];
	    }
	    return ...; // An array (in p) with the value of the integrands
	}
    };

    /**
     *   Class for evaluating the integrand associated with the
     *   Gamma-bar (advection like) coefficient in the Svensson
     *   transport. This coefficient only depends on the advection (A)
     *   and diffusion (D) coefficients, so there are two classes
     *   handeling the two input coefficients separately.
     */
    // YYY Work out a more descriptive name
    class SvenssonTransportGammaBarA : public SvenssonTransport<FVM::AdvectionTerm>{
	// Function for calculating the integrand associated to the
	// advection coefficient
	const reat_t* EvaluateIntegrand(len_t ir){
	    // The E field (array over nr) are taken from:
	    real_t *E = this->unknowns->GetUnknownData(this->EID);
	    // The rest of the quantities needed can be foun from:
	    // this->REFluid->Get...()
	    for( len_t i=0; i < this->np; i++){
		pBar = E[ir];
	    }
	    return ...; // An array (in p) with the value of the integrands
	}
    };
}

#include "SvensonTransport.tcc"

#endif/*_DREAM_SVENSSON_TRANSPORT_HPP*/


// XXX


