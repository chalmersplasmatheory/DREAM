#ifndef _DREAM_FVM_MATRIX_INVERTER_GMRES_HPP
#define _DREAM_FVM_MATRIX_INVERTER_GMRES_HPP

#include <petscksp.h>
#include "FVM/config.h"
#include "FVM/MatrixInverter.hpp"

namespace DREAM::FVM {
	class MIGMRES : public MatrixInverter {
    private:
        Vec x;

        PetscScalar *x_data;
        len_t xn;
	public:
		MIGMRES(const len_t);
        ~MIGMRES();

		virtual void Invert(Matrix*, Vec*, Vec*) override;

        void SetConvergenceTest(
            PetscErrorCode (*)(KSP, PetscInt, PetscReal, KSPConvergedReason*, void*),
            void*
        );
	};
}

#endif/*_DREAM_FVM_MATRIX_INVERTER_GMRES_HPP*/
