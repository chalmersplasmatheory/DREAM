#ifndef _DREAM_FVM_MATRIX_INVERTER_GMRES_HPP
#define _DREAM_FVM_MATRIX_INVERTER_GMRES_HPP

#include <petscksp.h>
#include <vector>
#include "FVM/config.h"
#include "FVM/MatrixInverter.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM {
	class MIGMRES : public MatrixInverter {
    private:
        Vec x = nullptr;

        PetscScalar *x_data = nullptr;
        len_t xn;

        PetscInt *blocks = nullptr;
        len_t nBlocks;
	public:
		MIGMRES(const len_t, std::vector<len_t>&, UnknownQuantityHandler*,
            PetscErrorCode (*)(KSP, PetscInt, PetscReal, KSPConvergedReason*, void*),
            void*
        );
        ~MIGMRES();

		virtual void Invert(Matrix*, Vec*, Vec*) override;
        void SetConvergenceTest(
            PetscErrorCode (*)(KSP, PetscInt, PetscReal, KSPConvergedReason*, void*),
            void*
        );
	};
}

#endif/*_DREAM_FVM_MATRIX_INVERTER_GMRES_HPP*/
