#ifndef _DREAM_FVM_MATRIX_INVERTER_MUMPS_HPP
#define _DREAM_FVM_MATRIX_INVERTER_MUMPS_HPP

#include <petscksp.h>
#include "FVM/config.h"
#include "FVM/MatrixInverter.hpp"

namespace DREAM::FVM {
	class MIMUMPS : public MatrixInverter {
    private:
        KSP ksp;
        Vec x;

        PetscScalar *x_data;
        len_t xn;
	public:
		MIMUMPS(const len_t);
        ~MIMUMPS();

		virtual void Invert(Matrix*, Vec*, Vec*) override;
	};
}

#endif/*_DREAM_FVM_MATRIX_INVERTER_MUMPS_HPP*/
