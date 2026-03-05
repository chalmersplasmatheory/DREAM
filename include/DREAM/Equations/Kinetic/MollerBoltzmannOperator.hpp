#ifndef _DREAM_EQUATIONS_MOLLER_BOLTZMANN_OPERATOR_HPP
#define _DREAM_EQUATIONS_MOLLER_BOLTZMANN_OPERATOR_HPP

#include "DREAM/Equations/Kinetic/MollerDeltaAngleKernel.hpp"
#include "DREAM/Equations/Kinetic/MollerEnergyKernel.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/config.h"

namespace DREAM {

/**
 * Generalized knock on Boltzmann source term, configured
 * for free-free large-angle Møller collisions.
 *
 * The term depends on unknowns:
 *  - n_tot: time-implicit
 *  - f_primary (f_hot or f_re): time explicit
 */
class MollerBoltzmannOperator : public FVM::EquationTerm {
   private:
    const FVM::Grid *gridPrimary = nullptr;
    FVM::UnknownQuantityHandler *unknowns = nullptr;

    len_t id_f_primary = 0;

    real_t scaleFactor = 1.0;

    // The two “helper” kernel owners.
    const MollerEnergyKernel *energyKernel;
    const MollerDeltaAngleKernel *angleKernel;

    // Assembled source vector on knock-on grid.
    real_t *sourceVector = nullptr;
    real_t t_source_rebuilt;

    // Scratch buffers reused in SetSourceVector().
    real_t *primaryWeights = nullptr;  // size: Np1P*NxiP, layout [k*NxiP + l]
    real_t *Cj = nullptr;  // size: NxiK

    void ValidateInputParameters() const;
    void ValidateGridAssumptions() const;

    void AllocateScratchBuffers();
    void AllocateSourceVector();
    void Deallocate();

    void BuildPrimaryWeights(len_t ir, const real_t *f_primary_ir, real_t *W_k_l) const;
    void SetSourceVector(const real_t *f_primary);

   public:
    MollerBoltzmannOperator(
        FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, FVM::UnknownQuantityHandler *,
        len_t id_f_primary, const MollerEnergyKernel *, const MollerDeltaAngleKernel *,
        real_t scaleFactor = 1.0
    );
    ~MollerBoltzmannOperator();

    virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *) override;
    virtual bool GridRebuilt() override;
    virtual void SetVectorElements(real_t *, const real_t *) override;
    virtual bool SetJacobianBlock(
        const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *x
    ) override;
    virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
};

}  // namespace DREAM

#endif /*_DREAM_EQUATIONS_MOLLER_BOLTZMANN_OPERATOR_HPP*/
