#ifndef _DREAM_EQUATIONS_FLUID_NBI_ION_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_NBI_ION_TERM_HPP

//#include "FVM/Equation/IonEquationTerm.hpp"
#include "DREAM/NBIHandler.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"  
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/RadialGrid.hpp"
#include "DREAM/IonHandler.hpp"





namespace DREAM {

class NBIIonTerm : public IonEquationTerm<FVM::EquationTerm> {

    private:
    len_t id_ncold, id_Tcold, id_ion_density, id_ion_temperature, id_ni;
    len_t iz;
    NBIHandler *handler;  
    IonHandler *ions;
    FVM::RadialGrid *radialGrid;
    len_t nr;
    FVM::UnknownQuantityHandler *unknowns;
public:
    
    NBIIonTerm(NBIHandler* h, FVM::Grid* grid, IonHandler* ionHandler,
               FVM::UnknownQuantityHandler* unknowns, const len_t iIon);
    
    virtual ~NBIIonTerm() override = default;

    virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler *unknowns) override; 
    virtual bool SetCSJacobianBlock(const len_t uqtyId, const len_t derivId,
                        FVM::Matrix *jac, const real_t *x,
                        const len_t iIon, const len_t Z0, const len_t rOffset
    ) override;
    virtual void SetCSMatrixElements(
        FVM::Matrix *mat, real_t *rhs,
        const len_t iIon, const len_t Z0, const len_t rOffset
    ) override;
    virtual void SetCSVectorElements(
        real_t *vec, const real_t *x,
        const len_t iIon, const len_t Z0, const len_t rOffset
    ) override;
    // Required by EquationTerm (pure in base)
    len_t GetNumberOfNonZerosPerRow() const override { return 30; }
    len_t GetNumberOfNonZerosPerRow_jac() const override { return 30; }


};


} 
#endif /*_DREAM_EQUATIONS_FLUID_NBI_ION_TERM_HPP */



