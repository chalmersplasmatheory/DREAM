#ifndef _DREAM_FVM_MATRIX_INVERTER_KSP_HPP
#define _DREAM_FVM_MATRIX_INVERTER_KSP_HPP

#include <petscksp.h>
#include "FVM/config.h"
#include "FVM/MatrixInverter.hpp"

namespace DREAM::FVM {
	class MIKSP : public MatrixInverter {
    private:
        KSP ksp;
        Vec x;

        PetscScalar *x_data;
        len_t xn;
	public:
		MIKSP(const len_t);
        ~MIKSP();

		virtual void Invert(Matrix*, Vec*, Vec*) override;

        void SetConvergenceTest(
            PetscErrorCode (*)(KSP, PetscInt, PetscReal, KSPConvergedReason*, void*),
            void*
        );
	};
}

#endif/*_DREAM_FVM_MATRIX_INVERTER_KSP_HPP*/
