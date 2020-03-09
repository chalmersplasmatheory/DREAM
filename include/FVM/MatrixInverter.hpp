#ifndef _TQS_FVM_MATRIX_INVERTER_HPP
#define _TQS_FVM_MATRIX_INVERTER_HPP

#include <petscvec.h>
#include "FVM/Matrix.hpp"

namespace TQS::FVM {
	class MatrixInverter {
	protected:
		Vec *solution = nullptr;
	public:
		MatrixInverter() {}

		virtual void Invert(Matrix*, Vec*, Vec*) = 0;
	};
}

#endif/*_FVM_MATRIX_INVERTER_HPP*/
