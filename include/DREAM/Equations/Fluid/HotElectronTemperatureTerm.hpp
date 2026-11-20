#ifndef _DREAM_EQUATION_HOT_ELECTRON_TEMPERATURE_TERM_HPP
#define _DREAM_EQUATION_HOT_ELECTRON_TEMPERATURE_TERM_HPP

namespace DREAM {
	class HotElectronHeatTerm : public EvaluableEquationTerm {
	protected:
		FVM::UnknownQuantityHandler *unknowns;

	public:
		HotElectronHeatTerm(FVM::Grid*, FVM::UnknownQuantityHandler*);

		virtual len_t GetNumberOfNonZerosPerRow() const override { return 2; }
		virtual len_t GetNumberOfNonZerosPerRow_jac() const override
		{ return 3; }

		virtual void EvaluableTransform(real_t*) override;
		virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
		virtual void SetVectorElements(real_t*, const real_t*) override;

		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
		virtual bool SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *x) override
		{ return this->DiagonalTerm::SetJacobianBlock(uqtyId, derivId, jac, x); }
	};
}

#endif/*_DREAM_EQUATION_HOT_ELECTRON_TEMPERATURE_TERM_HPP*/
