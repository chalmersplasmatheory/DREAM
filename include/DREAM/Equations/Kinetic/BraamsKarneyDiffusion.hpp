#ifndef _DREAM_EQUATIONS_KINETIC_BRAAMS_KARNEY_DIFFUSION_HPP
#define _DREAM_EQUATIONS_KINETIC_BRAAMS_KARNEY_DIFFUSION_HPP

#include "FVM/config.h"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"

namespace DREAM {
    class BraamsKarneyDiffusion
        : public FVM::DiffusionTerm {
    private:
		struct DCoeffs {
			real_t **dd11=nullptr, **dd12=nullptr,
				**dd21=nullptr, **dd22=nullptr;

			void init(len_t nr, len_t *n1, len_t *n2);
			~DCoeffs();
		};

		static constexpr int dds_n1 = 10, dds_n2 = 6;
		static constexpr int n1_offsets[dds_n1] = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4};
		static constexpr int n2_offsets[dds_n2] = {-3, -2, -1, 0, 1, 2};
		DCoeffs dds[dds_n2][dds_n1];
        
        CoulombLogarithm *lnLambda;
        len_t id_upsilon_1, id_upsilon_2;
        virtual void SetPartialDiffusionTerm(len_t derivId, len_t nMultiples) override;
		virtual void SetPartialJacobianContribution(int_t diagonalOffset, jacobian_interp_mode set_mode, len_t n, FVM::Matrix *jac, const real_t *x, bool momentumGrid) override;

        template<typename T1, typename T2>
		void SetCoefficients(T1 psi, T2 phi,
							 real_t **d11, real_t **d12,
							 real_t **d21, real_t **d22);
    public:

		BraamsKarneyDiffusion(FVM::Grid *g, FVM::UnknownQuantityHandler *unknowns, CoulombLogarithm *lnLambda, len_t id_upsilon_1, len_t id_upsilon_2);
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_BRAAMS_KARNEY_DIFFUSION_HPP*/

