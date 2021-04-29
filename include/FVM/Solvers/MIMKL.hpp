#ifndef _DREAM_FVM_MATRIX_INVERTER_MKL_HPP
#define _DREAM_FVM_MATRIX_INVERTER_MKL_HPP

#include <petscksp.h>
#include "FVM/config.h"
#include "FVM/MatrixInverter.hpp"

namespace DREAM::FVM {
	class MIMKL : public MatrixInverter {
    private:
        bool verbose = false;
        Vec x;

        PetscScalar *x_data;
        len_t xn;
	public:
		MIMKL(const len_t, bool verbose=false);
        ~MIMKL();

		virtual void Invert(Matrix*, Vec*, Vec*) override;
	};
}

#endif/*_DREAM_FVM_MATRIX_INVERTER_MKL_HPP*/
