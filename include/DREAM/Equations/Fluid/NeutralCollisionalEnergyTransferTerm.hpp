#ifndef DREAM_NEUTRAL_COLLISIONAL_ENERGY_TRANSFER_TERM_HPP
#define DREAM_NEUTRAL_COLLISIONAL_ENERGY_TRANSFER_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include <string>

namespace DREAM {

class NeutralCollisionalEnergyTransferTerm : public FVM::EquationTerm {
private:
    FVM::UnknownQuantityHandler *unknowns;
    IonHandler *ions;
    len_t iz;
    len_t id_ions;
    len_t id_Wn;
    real_t mi;


    real_t getMomentumScatteringRadius(len_t species);
    void GetNeutralParameters(len_t ir, real_t&, real_t&, real_t&);
    real_t CalculateNeutralScatteringRate(len_t, len_t, len_t, real_t, real_t, real_t, real_t, real_t);


public:
    NeutralCollisionalEnergyTransferTerm(
        FVM::Grid*, len_t, FVM::UnknownQuantityHandler*, IonHandler*
    );
    virtual ~NeutralCollisionalEnergyTransferTerm() {}
    virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}
    virtual void SetVectorElements(real_t*, const real_t*) override;
    virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
    virtual len_t GetNumberOfNonZerosPerRow() const override {
      return 0; //TODO
  }
};


}

#endif
