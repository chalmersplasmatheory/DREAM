#ifndef _DREAM_EQUATIONS_KINETIC_BRAAMS_KARNEY_ADVECTION_HPP
#define _DREAM_EQUATIONS_KINETIC_BRAAMS_KARNEY_ADVECTION_HPP


#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "FVM/config.h"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class BraamsKarneyAdvection
        : public FVM::AdvectionTerm {
    private:
		struct DCoeffs {
			real_t **df1=nullptr, **df2=nullptr;

			void init(len_t nr, len_t *n1, len_t *n2);
			~DCoeffs();
		};

		static constexpr int dfs_n1 = 4, dfs_n2 = 2;
		static constexpr int n1_offsets[dfs_n1] = {-2, -1, 0, 1};
		static constexpr int n2_offsets[dfs_n2] = {-1, 0};
		DCoeffs dfs[dfs_n2][dfs_n1];

		template<typename T1>
		void SetCoefficients(T1 psi, bool overwrite, real_t **df1, real_t **df2);

        CoulombLogarithm *lnLambda;
        len_t id_pi_0, id_pi_1;

        virtual void SetPartialAdvectionTerm(len_t derivId, len_t nMultiples) override;
		virtual void SetPartialJacobianContribution(int_t diagonalOffset, jacobian_interp_mode set_mode, len_t n, FVM::Matrix *jac, const real_t *x, bool momentumGrid) override;

    public:
        BraamsKarneyAdvection(FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns, CoulombLogarithm *lnLambda, len_t id_pi_0, len_t id_pi_1);
        
        
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_BRAAMS_KARNEY_ADVECTION_HPP*/
