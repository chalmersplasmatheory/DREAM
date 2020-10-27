#ifndef _DREAM_FVM_MATRIX_INVERTER_MUMPS_HPP
#define _DREAM_FVM_MATRIX_INVERTER_MUMPS_HPP

#include <petscksp.h>
#include "FVM/config.h"
#include "FVM/MatrixInverter.hpp"

namespace DREAM::FVM {
	class MIMUMPS : public MatrixInverter {
    public:
        enum ICNTL {
            ICNTL_OUTPUT_STREAM_ERROR=1,
            ICNTL_OUTPUT_STREAM_DIAGNOSTIC=2,
            ICNTL_OUTPUT_STREAM_GLOBAL=3,
            ICNTL_PRINTING_LEVEL=4,

            // ...
            // TODO
            // ...

            ICNTL_OPENMP_THREADS=16,
            ICNTL_DETECT_NULL_PIVOT_ROWS=24
        };
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
