#ifndef _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_NBI_ELECTRON_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/NBIHandler.hpp"
#include "DREAM/IonHandler.hpp"

//#include "DREAM/UnknownQuantityHandler.hpp"


namespace DREAM {

class NBIElectronTerm : public FVM::EquationTerm {

private:
    NBIHandler *handler;    
    FVM::RadialGrid *radialGrid;
    len_t nr;
    FVM::UnknownQuantityHandler *unknowns;
    len_t id_ncold, id_Tcold, id_ion_density, id_ion_temperature ,id_ni;
    IonHandler* ions;
    const real_t* Qe_1;
    const real_t* Qe_2;
    const real_t* Qe_3;
public:
    NBIElectronTerm(NBIHandler *h, FVM::Grid *grid, FVM::UnknownQuantityHandler *unknowns, IonHandler* ionHandler);

    ~NBIElectronTerm() override = default;

    virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns) override;
    virtual void SetVectorElements(real_t *rhs, const real_t *x) override;
    virtual void SetMatrixElements(FVM::Matrix *mat, real_t *rhs) override;
    virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
    virtual len_t GetNumberOfNonZerosPerRow() const override { return 30; }
    virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 30; }


};


} 
#endif /*_DREAM_EQUATIONS_FLUID_NBI_ELECTRON_TERM_HPP */



