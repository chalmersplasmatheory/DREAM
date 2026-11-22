#ifndef _DREAM_ROSENBLUTH_PUTVINSKI_BOLTZMANN_HPP
#define _DREAM_ROSENBLUTH_PUTVINSKI_BOLTZMANN_HPP

#include <array>
#include <vector>
#include "FVM/Equation/EquationTerm.hpp"

namespace DREAM {

/**
 * Infinite-energy RP spectrum and previous-step off-grid primary density.
 * Momentum is in units of m_e c.
 */
class RosenbluthPutvinskiBoltzmann {
public:
    enum Region { HOTTAIL, RUNAWAY, RUNAWAY_DENSITY, ALL, ABOVE_GRID };

    RosenbluthPutvinskiBoltzmann(
        const FVM::Grid *hottail, const FVM::Grid *runaway,
        FVM::UnknownQuantityHandler *unknowns, real_t pCutoff, len_t poloidalOrder
    );

    void RebuildPrimaryDensities(FVM::UnknownQuantityHandler *unknowns);
    void GridRebuilt();

    const real_t *GetSource(Region region) const;
    real_t GetProduction(Region region) const { return production[region]; }
    const real_t *GetExternalDensity() const { return external.data(); }
    const real_t *GetRawExternalDensity() const { return externalRaw.data(); }

private:
    const FVM::Grid *hottail;
    const FVM::Grid *runaway;
    FVM::UnknownQuantityHandler *unknowns;
    len_t idNRe, idFRe = 0;
    real_t pCutoff;
    len_t poloidalOrder;
    std::vector<real_t> external, externalRaw;
    std::array<std::vector<real_t>, 2> source;
    std::array<real_t, 5> production{};
    real_t previousTime;
    bool densitiesValid = false;

    void BuildSource(Region, const FVM::Grid *);
};

/**
 * Infinite-energy source multiplying n_tot. Either the kinetic shape or its
 * fluid moment over the selected region. Both use the same cached primaries
 * and spectrum; diagnostics can evaluate each parent population separately.
 */
class RosenbluthPutvinskiBoltzmannTerm : public FVM::EquationTerm {
public:
    RosenbluthPutvinskiBoltzmannTerm(
        FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns,
        RosenbluthPutvinskiBoltzmann *source, RosenbluthPutvinskiBoltzmann::Region region,
        real_t scaleFactor
    );

    void Rebuild(real_t, real_t, FVM::UnknownQuantityHandler *) override;
    void SetVectorElements(real_t *vec, const real_t *nTot) override;
    void SetMatrixElements(FVM::Matrix *, real_t *) override;
    bool SetJacobianBlock(len_t uqtyId, len_t derivId, FVM::Matrix *, const real_t *nTot) override;
    len_t GetNumberOfNonZerosPerRow() const override { return 1; }

    void SetExternalSource(real_t *vec, const real_t *nTot) const;
    void GetExternalDensity(real_t *density) const;
    void GetRawExternalDensity(real_t *density) const;
    void GetAboveGridProduction(real_t *rate, const real_t *nTot) const;

private:
    FVM::UnknownQuantityHandler *unknowns;
    RosenbluthPutvinskiBoltzmann *source;
    RosenbluthPutvinskiBoltzmann::Region region;
    real_t scaleFactor;
    bool kinetic;
    template<typename Add>
    void Apply(const real_t *nTot, Add add) const;
};

} // namespace DREAM

#endif
