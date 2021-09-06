#ifndef _DREAM_FVM_MATRIX_INVERTER_GMRES_HPP
#define _DREAM_FVM_MATRIX_INVERTER_GMRES_HPP

#include <petscksp.h>
#include <vector>
#include "FVM/config.h"
#include "FVM/MatrixInverter.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM::FVM {
	class MIGMRES : public MatrixInverter {
    private:
        Vec x;

        PetscScalar *x_data;
        len_t xn;

        PetscInt *blocks;
        len_t nBlocks;
	public:
		MIGMRES(const len_t, std::vector<len_t>&, UnknownQuantityHandler*);
        ~MIGMRES();

		virtual void Invert(Matrix*, Vec*, Vec*) override;
	};
}

#endif/*_DREAM_FVM_MATRIX_INVERTER_GMRES_HPP*/
