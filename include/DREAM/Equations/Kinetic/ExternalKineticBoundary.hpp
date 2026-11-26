#ifndef _DREAM_EXTERNAL_KINETIC_BOUNDARY_HPP
#define _DREAM_EXTERNAL_KINETIC_BOUNDARY_HPP

#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>
#include "FVM/Equation/BoundaryConditions/PXiAdvectionDiffusionBoundaryCondition.hpp"

namespace DREAM {

/**
 * Two-way outer p boundary. Passing flux is assigned by its pitch sign
 * and trapped flux is shared between the two external inventories.
 * The kinetic equation sees the sum of the two transfers.
 */
class ExternalKineticBoundary : public FVM::BC::PXiAdvectionDiffusionBoundaryCondition {
public:
    enum Side { KINETIC, DENSITY, NEGATIVE_DENSITY, POSITIVE_DENSITY };
    enum Pole { POSITIVE, NEGATIVE, BOTH };

    ExternalKineticBoundary(
        FVM::Grid *equationGrid, FVM::Grid *distributionGrid,
        const FVM::Operator *, len_t idF, std::array<len_t, 2> idExternal,
        FVM::UnknownQuantityHandler *, Side side=KINETIC
    );

    bool GridRebuilt() override;
    bool Rebuild(real_t t, FVM::UnknownQuantityHandler *) override;
    len_t GetNumberOfNonZerosPerRow() const override;
    len_t GetNumberOfNonZerosPerRow_jac() const override;
    void AddToVectorElements(real_t *, const real_t *) override;
    bool AddToJacobianBlock(len_t, len_t, FVM::Matrix *, const real_t *) override;
    void AddToMatrixElements(FVM::Matrix *, real_t *) override;
    void AddToVectorElements_c(
        real_t *, const real_t *, const real_t *const*, const real_t *const*,
        const real_t *const*, const real_t *const*, const real_t *const*,
        const real_t *const*, const real_t *const*, const real_t *const*,
        jacobian_interp_mode=NO_JACOBIAN
    ) override;

    bool SetJacobianBlock(len_t, len_t, FVM::Matrix *, const real_t *) override { return false; }
    void SetMatrixElements(FVM::Matrix *, real_t *) override {}
    void SetVectorElements(real_t *, const real_t *) override {}

    // Physical rates, positive out of the kinetic grid. No minus sign
    // from the density equation is included in these diagnostics.
    void GetFlux(real_t *net, real_t *outward=nullptr, real_t *inward=nullptr, Pole=BOTH) const;
    void GetRequestedInflow(real_t *, Pole=BOTH) const;

private:
    using PoleValues = std::array<real_t, 2>;
    struct Trace {
        real_t value, slope;
        // Derivatives with respect to the penultimate and last cells.
        std::array<real_t, 2> dValue, dSlope;
    };
    struct Face {
        real_t flux=0, weight=0, kineticScale=0;
        PoleValues fraction{};
        real_t advectionValue=0, slope=0, pitchSlope=0;
        // Incoming trace comes from the last accepted distribution, evaluated
        // with the transport coefficients of the current Newton iterate.
        real_t acceptedFlux=0, inflow=0;
        real_t inflowAdvectionValue=0, inflowSlope=0, inflowPitchSlope=0;
        // Columns 2*j and 2*j+1 are the last two momentum cells at pitch j.
        std::vector<std::pair<len_t, real_t>> derivative;
    };
    struct Reservoir {
        real_t scale=0, dDensity=0, dRate=0;
    };
    using Reservoirs = std::array<Reservoir, 2>;

    FVM::Grid *distributionGrid;
    len_t idF;
    std::array<len_t, 2> idExternal;
    FVM::UnknownQuantityHandler *unknowns;
    Side side;
    real_t dt=0;
    // Retain the two accepted boundary cells through post-step output saving.
    std::vector<std::vector<real_t>> acceptedBoundaryValues;
    std::vector<real_t> positiveFraction;
    std::vector<len_t> pitchOffset;

    void BuildPoleFractions();
    Trace Reconstruct(real_t previous, real_t last, real_t distance, real_t halfWidth) const;
    std::vector<Face> EvaluateFaces(len_t ir, len_t offset, const real_t *f) const;
    std::vector<Face> EvaluateRawFaces(len_t ir, len_t offset, const real_t *f,
        bool accepted=false) const;
    static PoleValues RequestedInflow(const std::vector<Face>&);
    Reservoirs LimitInflow(len_t ir, const PoleValues&) const;
    Reservoir LimitInflow(real_t density, real_t rate) const;
    real_t AppliedFlux(const Face&, const Reservoirs&) const;
    void AddDistributionJacobian(len_t ir, len_t offset, const std::vector<Face>&, FVM::Matrix*) const;
    bool IncludesPole(len_t pole) const {
        return (side != POSITIVE_DENSITY || pole == POSITIVE)
            && (side != NEGATIVE_DENSITY || pole == NEGATIVE);
    }
    len_t Row(len_t ir, len_t offset, len_t j) const;
    real_t Scale(const Face&) const;
};

}
#endif
