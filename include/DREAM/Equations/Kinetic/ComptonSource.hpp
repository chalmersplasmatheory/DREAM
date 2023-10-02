#ifndef _DREAM_EQUATIONS_COMPTON_SOURCE_HPP
#define _DREAM_EQUATIONS_COMPTON_SOURCE_HPP

#include "DREAM/Equations/FluidSourceTerm.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include <limits>

namespace DREAM {
    class ComptonSource
        : public FluidSourceTerm {
    
    public:
        enum SourceMode{
            SOURCE_MODE_FLUID,
            SOURCE_MODE_KINETIC
        };
    private:
        len_t id_ne, id_Eterm, id_ntot;
        
        real_t *source;
        FVM::Interpolator1D *comptonPhotonFlux;
        real_t photonFlux, pLower;
        real_t pc, scaleFactor;
        
        len_t limit;
        gsl_integration_workspace * wp;
        gsl_integration_workspace * wpOut;
        
        SourceMode sourceMode;
        
    	RunawayFluid *REFluid;
    protected:
        static real_t integrand(real_t, void * );
        static real_t innerIntegrand(real_t, void * );
        static real_t fluidIntegrand(real_t, void * );
        static real_t integratedPhotonEnergySpectrum(real_t , void * );
        
        virtual real_t GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId) override;
    public:
        struct intparams{
            len_t limit; 
            gsl_integration_workspace * wp;
        };
        ComptonSource(FVM::Grid*, FVM::UnknownQuantityHandler*, FVM::Interpolator1D*, real_t=0., real_t=0., SourceMode sm = SOURCE_MODE_KINETIC, RunawayFluid* REFluid=nullptr);
        
        virtual real_t GetSourceFunction(len_t ir, len_t i, len_t j) override;

        real_t EvaluateSource(len_t ir, len_t i, len_t j);
        static real_t EvaluateTotalComptonNumber(real_t pLower, intparams * params, intparams * paramsOut, real_t pUpper=std::numeric_limits<real_t>::infinity());
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_COMPTON_SOURCE_HPP*/
