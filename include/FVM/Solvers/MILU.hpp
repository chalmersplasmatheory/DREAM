#ifndef _TQS_FVM_MATRIX_INVERTER_LU_HPP
#define _TQS_FVM_MATRIX_INVERTER_LU_HPP

#include <petscksp.h>
#include "FVM/config.h"
#include "FVM/MatrixInverter.hpp"

namespace TQS::FVM {
	class MILU : public MatrixInverter {
    private:
        KSP ksp;
        Vec x;

        PetscScalar *x_data;
        len_t xn;
	public:
		MILU(const len_t);
        ~MILU();

		virtual void Invert(Matrix*, Vec*, Vec*) override;
	};
}

#endif/*_TQS_FVM_MATRIX_INVERTER_LU_HPP*/
