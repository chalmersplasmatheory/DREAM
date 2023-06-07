#ifndef _DREAM_ION_SOURCE_BOUNDARY_CONDITION_HPP
#define _DREAM_ION_SOURCE_BOUNDARY_CONDITION_HPP

#include "DREAM/MultiInterpolator1D.hpp"
#include "DREAM/Equations/Fluid/IonBoundaryCondition.hpp"
#include "FVM/Equation/BoundaryCondition.hpp"

namespace DREAM {
	class IonSourceBoundaryCondition : public IonBoundaryCondition<FVM::BC::BoundaryCondition> {
	protected:
		MultiInterpolator1D *dndt = nullptr;
		real_t time = 0;

	public:
		IonSourceBoundaryCondition(
			FVM::Grid*, IonHandler*, MultiInterpolator1D*,
			const len_t iIon
		);

		virtual ~IonSourceBoundaryCondition();

		virtual bool Rebuild(const real_t, FVM::UnknownQuantityHandler*) override;

		/* Derivative of a constant is zero, so the jacobian
		 * for this term vanishes */
		virtual bool AddCSToJacobianBlock(
			const len_t, const len_t, FVM::Matrix*, const real_t*,
			const len_t, const len_t, const len_t
		) override {return false;}
		virtual bool SetCSJacobianBlock(
			const len_t, const len_t, FVM::Matrix*, const real_t*,
			const len_t, const len_t, const len_t
		) override {return false;}

		virtual void AddCSToMatrixElements(
			FVM::Matrix*, real_t*, const len_t iIon,
			const len_t Z0, const len_t rOffset
		) override;
		virtual void SetCSMatrixElements(
			FVM::Matrix*, real_t*, const len_t,
			const len_t, const len_t
		) override {}	// No overwriting of elements done

		virtual void AddCSToVectorElements(
			real_t*, const real_t*, const len_t iIon,
			const len_t Z0, const len_t rOffset
		) override;
		virtual void SetCSVectorElements(
			real_t*, const real_t*, const len_t,
			const len_t, const len_t
		) override {}	// No overwriting of elements done
	};
}

#endif/*_DREAM_ION_SOURCE_BOUNDARY_CONDITION_HPP*/
