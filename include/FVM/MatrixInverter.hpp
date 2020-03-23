#ifndef _DREAM_FVM_MATRIX_INVERTER_HPP
#define _DREAM_FVM_MATRIX_INVERTER_HPP

#include <petscvec.h>
#include "FVM/Matrix.hpp"

namespace DREAM::FVM {
	class MatrixInverter {
	protected:
		Vec *solution = nullptr;
	public:
		MatrixInverter() {}

		virtual void Invert(Matrix*, Vec*, Vec*) = 0;
	};
}

#endif/*_FVM_MATRIX_INVERTER_HPP*/
