#ifndef _DREAM_DIAGONAL_PRECONDITIONER_HPP
#define _DREAM_DIAGONAL_PRECONDITIONER_HPP

#include <unordered_map>
#include <vector>
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class DiagonalPreconditioner {
	public:
		static const real_t
			DENSITY_SCALE,
			CURRENT_SCALE,
			ENERGY_SCALE,
			FLUX_SCALE,
			TEMPERATURE_SCALE,
			RUNAWAY_FRACTION;

    private:
        FVM::UnknownQuantityHandler *uqh;
        std::vector<len_t> nontrivials;

        // Vectors representing diagonal scaling matrices and their inverses
        Vec iuqn, eqn;

        std::unordered_map<len_t, real_t> uqn_scales;   // Scaling factors for unknowns
        std::unordered_map<len_t, real_t> eqn_scales;   // Scaling factors for equations

    public:
        DiagonalPreconditioner(
            FVM::UnknownQuantityHandler*, const std::vector<len_t>&
        );
        ~DiagonalPreconditioner();

		real_t GetEquationScale(const len_t i) { return this->eqn_scales[i]; }
		real_t GetUnknownScale(const len_t i) { return this->uqn_scales[i]; }

        void Build();
        void SetEquationScale(const len_t, const real_t);
        void SetUnknownScale(const len_t, const real_t);
        void SetDefaultScalings();

        void RescaleMatrix(FVM::Matrix*);
        void RescaleRHSVector(Vec);

        void UnscaleUnknownVector(Vec);
    };
}

#endif/*_DREAM_DIAGONAL_PRECONDITIONER_HPP*/
