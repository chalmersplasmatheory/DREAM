/**
 * Implementation of the p*nu_s friction term in the kinetic equation.
 */

#include "DREAM/Equations/Kinetic/SlowingDownTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
SlowingDownTerm::SlowingDownTerm(
    FVM::Grid *g, CollisionQuantityHandler *cqh,
    enum OptionConstants::momentumgrid_type mgtype,
    FVM::UnknownQuantityHandler *unknowns,
    bool withKineticIonJacobian
) : FVM::AdvectionTerm(g) {

    this->gridtype = mgtype;
    this->nuS = cqh->GetNuS();
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD));
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
    if(withKineticIonJacobian)
        AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES));
}

/**
 * Build the coefficients of this advection term.
 */
void SlowingDownTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*){
    real_t *const* nu_s_f1 = nuS->GetValue_f1();
    real_t *const* nu_s_f2 = nuS->GetValue_f2();
  
    bool gridtypePXI        = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
    bool gridtypePPARPPERP  = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = grid->GetMomentumGrid(ir);
                
        if (gridtypePXI || gridtypePPARPPERP) 
            for (len_t j = 0; j < n2[ir]; j++) 
                for (len_t i = 0; i < n1[ir]+1; i++) 
                    F1(ir, i, j) -= mg->GetP1_f(i) * nu_s_f1[ir][j*(n1[ir]+1)+i];
                
        if (gridtypePPARPPERP) 
            for (len_t j = 0; j < n2[ir]+1; j++) 
                for (len_t i = 0; i < n1[ir]; i++) 
                    F2(ir, i, j) -= mg->GetP2_f(j) * nu_s_f2[ir][j*n1[ir]+i];
                
        if (gridtypePXI){    
            real_t p3nuSAtZero = nuS->GetP3NuSAtZero(ir);
            for(len_t j=0; j< n2[ir]; j++)
                F1PSqAtZero(ir,j) -= p3nuSAtZero;
        }
    }

}


// Set jacobian of the advection coefficients for this advection term
void SlowingDownTerm::SetPartialAdvectionTerm(len_t derivId, len_t nMultiples){
    ResetDifferentiationCoefficients();

    const len_t nr = grid->GetNr();
 
    // Get jacobian of the collision frequency
    const real_t *dNuS_f1 = nuS->GetUnknownPartialContribution(derivId, FVM::FLUXGRIDTYPE_P1);
    const real_t *dNuS_f2 = nuS->GetUnknownPartialContribution(derivId, FVM::FLUXGRIDTYPE_P2);
    real_t *dp3nuSAtZero  = nuS->GetPartialP3NuSAtZero(derivId); 
    bool gridtypePXI        = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
    bool gridtypePPARPPERP  = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);

    len_t offset1 = 0;
    len_t offset2 = 0;

    // Set partial advection coefficients
    for(len_t n=0; n<nMultiples; n++)
        for (len_t ir = 0; ir < nr; ir++) {
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = n1[ir];
        const len_t np2 = n2[ir];
            if (gridtypePXI || gridtypePPARPPERP) 
                for (len_t j = 0; j < np2; j++) 
                    for (len_t i = 0; i < np1+1; i++) 
                        dF1(ir,i,j,n) = -mg->GetP1_f(i) * dNuS_f1[offset1 + j*(np1+1) + i];

            if (gridtypePPARPPERP) 
                for (len_t j = 0; j < np2+1; j++) 
                    for (len_t i = 0; i < np1; i++) 
                        dF2(ir,i,j,n) = -mg->GetP2_f(j) * dNuS_f2[offset2 + j*np1 + i];

            if (gridtypePXI == (dp3nuSAtZero != nullptr)) {
                real_t dp3nuS = dp3nuSAtZero[nr*n+ir];
                for (len_t j=0; j<np2; j++)
                    dF1PSqAtZero(ir,j,n) = -dp3nuS;
            }
            offset1 += (np1+1)*np2;
            offset2 += np1*(np2+1);
        }
    delete [] dp3nuSAtZero;
}
