#include "DREAM/Equations/Fluid/IonNeutralAdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"

using namespace DREAM;

IonNeutralAdvectionTerm::IonNeutralAdvectionTerm(FVM::Grid *grid, IonHandler *ions, const len_t iIon,
    bool allocCoefficients, FVM::AdvectionInterpolationCoefficient::adv_interpolation intp,
    OptionConstants::adv_jacobian_mode jac_mode, len_t id, real_t damping_factor) : IonNeutralAdvectionDiffusionTerm<FVM::AdvectionTerm>(grid, ions, iIon, allocCoefficients){
        SetAdvectionInterpolationMethod(intp, jac_mode, id, damping_factor);
    }
    
IonNeutralAdvectionTerm::~IonNeutralAdvectionTerm(){}

void IonNeutralAdvectionTerm::Rebuild(const real_t t, const real_t dt , FVM::UnknownQuantityHandler* uqty){
    this->SetCoeffs(t);
    
    this->FVM::AdvectionTerm::RebuildFluxLimiterDamping(t, dt);
    this->FVM::AdvectionTerm::RebuildInterpolationCoefficients(uqty, nullptr, nullptr, nullptr);
}
		
		
		
