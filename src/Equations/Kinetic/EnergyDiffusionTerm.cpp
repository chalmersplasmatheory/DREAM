/**
 * Implementation of the nu_parallel collisional energy diffusion term in the kinetic equation.
 */

#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/Kinetic/EnergyDiffusionTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
EnergyDiffusionTerm::EnergyDiffusionTerm(FVM::Grid *g, CollisionQuantityHandler *cqh, 
    enum OptionConstants::momentumgrid_type mgtype, FVM::UnknownQuantityHandler *unknowns,
    bool withKineticIonJacobian)
    : FVM::DiffusionTerm(g) {
        this->gridtype = mgtype;
        this->nuPar    = cqh->GetNuPar();
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD));
    AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
    if(withKineticIonJacobian && !(nuPar->GetSettings()->screened_diffusion==OptionConstants::COLLQTY_SCREENED_DIFFUSION_MODE_ZERO))
        AddUnknownForJacobian(unknowns, unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES));
}

/**
 * Build the coefficients of this diffusion term.
 */
void EnergyDiffusionTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *){ 
    real_t *const* nu_par_f1 = nuPar->GetValue_f1();
    real_t *const* nu_par_f2 = nuPar->GetValue_f2();

    bool gridtypePXI       = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
    bool gridtypePPARPPERP = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        real_t xi0_f1, xi0_f2;        
        for (len_t j = 0; j < np2; j++)
            for (len_t i = 0; i < np1+1; i++) {
                if (gridtypePXI)
                    D11(ir, i, j) += nu_par_f1[ir][j*(np1+1)+i];
                else if (gridtypePPARPPERP){
                    xi0_f1 = mg->GetXi0_f1(i,j);
                    D11(ir, i, j) += xi0_f1*xi0_f1*nu_par_f1[ir][j*(np1+1)+i];
                    D12(ir, i,j ) += xi0_f1*sqrt(1-xi0_f1*xi0_f1)*nu_par_f1[ir][j*(np1+1)+i];
                } 
            }
        if (gridtypePPARPPERP) 
            for (len_t j = 0; j < np2+1; j++) 
                for (len_t i = 0; i < np1; i++) {
                    xi0_f2 = mg->GetXi0_f2(i,j);
                    D22(ir, i, j) += (1-xi0_f2*xi0_f2) * nu_par_f2[ir][j*np1+i];
                    D21(ir, i, j) += xi0_f2*sqrt(1-xi0_f2*xi0_f2) * nu_par_f2[ir][j*np1+i];
                }
    }
}




// Set jacobian of the diffusion coefficients for this diffusion term
void EnergyDiffusionTerm::SetPartialDiffusionTerm(len_t derivId, len_t nMultiples){
    ResetDifferentiationCoefficients();     
    const real_t *dNuPar_f1 = nuPar->GetUnknownPartialContribution(derivId, FVM::FLUXGRIDTYPE_P1);
    const real_t *dNuPar_f2 = nuPar->GetUnknownPartialContribution(derivId, FVM::FLUXGRIDTYPE_P2);

    bool gridtypePXI       = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
    bool gridtypePPARPPERP = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);

    len_t offset1 = 0;
    len_t offset2 = 0;
    for(len_t n = 0; n < nMultiples; n++)
        for (len_t ir = 0; ir < nr; ir++) {
            auto *mg = grid->GetMomentumGrid(ir);
            const len_t np1 = mg->GetNp1();
            const len_t np2 = mg->GetNp2();

            real_t xi0_f1, xi0_f2;
            for (len_t j = 0; j < np2; j++)
                for (len_t i = 0; i < np1+1; i++){
                    if (gridtypePXI)
                        dD11(ir, i, j, n) += dNuPar_f1[offset1 + j*(np1+1) + i];
                    else if (gridtypePPARPPERP){
                        xi0_f1 = mg->GetXi0_f1(i,j);
                        dD11(ir, i, j, n) += xi0_f1*xi0_f1 * dNuPar_f1[offset1 + j*(np1+1) + i];
                        dD12(ir, i,j, n)  += xi0_f1*sqrt(1-xi0_f1*xi0_f1) * dNuPar_f1[offset1 + j*(np1+1) + i];
                    }
                }

            if(gridtypePPARPPERP) 
                for (len_t j = 0; j < np2+1; j++) 
                    for (len_t i = 0; i < np1; i++) {
                        xi0_f2 = mg->GetXi0_f2(i,j);
                        dD22(ir, i, j, n) += (1-xi0_f2*xi0_f2) * dNuPar_f2[offset2 + j*np1 + i];
                        dD21(ir, i, j, n) += xi0_f2*sqrt(1-xi0_f2*xi0_f2) * dNuPar_f2[offset2 + j*np1 + i];
                    }
            offset1 += (np1+1)*np2;
            offset2 += np1*(np2+1);
        }
}


