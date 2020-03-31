#ifndef _DREAM_FVM_ADVECTION_DIFFUSION_TERM_HPP
#define _DREAM_FVM_ADVECTION_DIFFUSION_TERM_HPP

#include <vector>
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Equation/EquationTerm.hpp"


namespace DREAM::FVM {
    class AdvectionDiffusionTerm : public AdvectionTerm, public DiffusionTerm {
    public:
        enum advdiff_interpolation {
            AD_INTERP_CENTRED,
            AD_INTERP_BACKWARD,
            AD_INTERP_FORWARD
        };

    private:
        std::vector<AdvectionTerm*> advectionterms;
        std::vector<DiffusionTerm*> diffusionterms;

        enum advdiff_interpolation interpolationMethod = AD_INTERP_CENTRED;

    public:
        AdvectionDiffusionTerm(Grid *g, enum advdiff_interpolation intp)
            : AdvectionTerm(g, true), DiffusionTerm(g), interpolationMethod(intp) {}

        void Add(AdvectionTerm*);
        void Add(DiffusionTerm*);

        virtual void Rebuild(const real_t) override;
        void RebuildInterpolationCoefficients();
        void SetInterpolationCoefficientValues(const real_t);

    };
}

#endif/*_DREAM_FVM_ADVECTION_DIFFUSION_TERM_HPP*/
