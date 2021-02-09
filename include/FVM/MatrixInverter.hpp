#ifndef _DREAM_FVM_MATRIX_INVERTER_HPP
#define _DREAM_FVM_MATRIX_INVERTER_HPP

#include <petscksp.h>
#include <petscvec.h>
#include "FVM/Matrix.hpp"

namespace DREAM::FVM {
	class MatrixInverter {
	protected:
		Vec *solution = nullptr;
        KSP ksp;
	public:
		MatrixInverter() {}
        virtual ~MatrixInverter() {}

		virtual void Invert(Matrix*, Vec*, Vec*) = 0;

        virtual void PrintInfo();
	};
}

#endif/*_FVM_MATRIX_INVERTER_HPP*/
