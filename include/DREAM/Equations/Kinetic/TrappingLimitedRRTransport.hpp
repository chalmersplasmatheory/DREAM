#ifndef _DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP
#define _DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP

#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/Kinetic/RechesterRosenbluthTransport.hpp"

namespace DREAM {
class TrappingLimitedRRTransport :
public RechesterRosenbluthTransport {
protected:
FVM::UnknownQuantityHandler *unknowns;
PitchScatterFrequency *nuD;
CollisionQuantity::collqty_settings *collSettings;
real_t currentTime = 0;

len_t id_ncold, id_Tcold, id_ni;

real_t EvaluateNuDOnRadialFluxGrid(len_t, real_t) const;
real_t EvaluatePartialNuDOnRadialFluxGrid(len_t, real_t, len_t, len_t) const;

virtual void SetPartialDiffusionTerm(len_t, len_t) override;

public:
TrappingLimitedRRTransport(
FVM::Grid*, enum OptionConstants::momentumgrid_type,
FVM::Interpolator1D*, CollisionQuantityHandler*,
FVM::UnknownQuantityHandler*
);
virtual ~TrappingLimitedRRTransport();

virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
};
}

#endif/*_DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP*/
