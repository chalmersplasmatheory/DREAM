#ifndef _DREAM_FVM_MATRIX_INVERTER_LU_HPP
#define _DREAM_FVM_MATRIX_INVERTER_LU_HPP

#include <petscksp.h>
#include "FVM/config.h"
#include "FVM/MatrixInverter.hpp"

namespace DREAM::FVM {
	class MILU : public MatrixInverter {
    private:
        len_t xn;
	public:
		MILU(const len_t);
        ~MILU();

		virtual void Invert(Matrix*, Vec*, Vec*) override;
	};
}

#endif/*_DREAM_FVM_MATRIX_INVERTER_LU_HPP*/
