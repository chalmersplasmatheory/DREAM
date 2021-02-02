#ifndef _DREAM_FVM_MATRIX_INVERTER_SUPERLU_HPP
#define _DREAM_FVM_MATRIX_INVERTER_SUPERLU_HPP

#include <petscksp.h>
#include "FVM/config.h"
#include "FVM/MatrixInverter.hpp"

namespace DREAM::FVM {
	class MISuperLU : public MatrixInverter {
    private:
        Vec x;

        PetscScalar *x_data;
        len_t xn;
	public:
		MISuperLU(const len_t);
        ~MISuperLU();

		virtual void Invert(Matrix*, Vec*, Vec*) override;
	};
}

#endif/*_DREAM_FVM_MATRIX_INVERTER_SUPERLU_HPP*/
