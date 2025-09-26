#ifndef _DREAM_EQUATIONS_FLUID_NBI_ION_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_NBI_ION_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/NBIHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "DREAM/IonHandler.hpp"





namespace DREAM {

class NBIIonTerm : public FVM::EquationTerm {

    private:
    len_t id_ncold, id_Tcold, id_ion_density, id_ion_temperature;
    len_t iz;
    NBIHandler *handler;  
    IonHandler *ions;
    FVM::RadialGrid *radialGrid;
    len_t nr;
    FVM::UnknownQuantityHandler *unknowns;
public:
    
    NBIIonTerm(NBIHandler *h, FVM::Grid *grid, IonHandler *ions, FVM::UnknownQuantityHandler *unknowns, len_t iz);
    
    ~NBIIonTerm() override = default;

    virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns) override; 
    void SetVectorElements(real_t *rhs, const real_t *x) override;
    void SetMatrixElements(FVM::Matrix *mat, real_t *rhs) override;
    bool SetJacobianBlock(const len_t uqtyId, const len_t derivId,FVM::Matrix *jac, const real_t *x) override;
    len_t GetNumberOfNonZerosPerRow() const override { return 30; }
    len_t GetNumberOfNonZerosPerRow_jac() const override { return 30; }

};


} 
#endif /*_DREAM_EQUATIONS_FLUID_NBI_ION_TERM_HPP */



