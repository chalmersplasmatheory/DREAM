/**
 * Implementation of the p*nu_s friction term in the kinetic equation.
 */

#include <softlib/SFile.h>
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/Kinetic/SlowingDownTerm.hpp"


using namespace DREAM;

/**
 * Constructor.
 */
SlowingDownTerm::SlowingDownTerm(
    FVM::Grid *g, CollisionQuantityHandler *cqh,
    enum OptionConstants::momentumgrid_type mgtype
) : FVM::AdvectionTerm(g) {

    this->gridtype = mgtype;
    this->nuS = cqh->GetNuS();
}

/**
 * Build the coefficients of this advection term.
 */
void SlowingDownTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler *unknowns){
    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    nzs      = unknowns->GetUnknown(id_ni)->NumberOfMultiples();

    real_t *const* nu_s_f1 = nuS->GetValue_f1();
    real_t *const* nu_s_f2 = nuS->GetValue_f2();
  
    bool gridtypePXI        = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
    bool gridtypePPARPPERP  = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);

    for (len_t ir = 0; ir < nr; ir++) {
        auto *mg = grid->GetMomentumGrid(ir);
                
        if (gridtypePXI || gridtypePPARPPERP) {
            for (len_t j = 0; j < n2[ir]; j++) {
                for (len_t i = 0; i < n1[ir]+1; i++) {
                    F1(ir, i, j) -= mg->GetP1_f(i) * nu_s_f1[ir][j*(n1[ir]+1)+i];
                }
            }
        }

        if (gridtypePPARPPERP) {
            for (len_t j = 0; j < n2[ir]+1; j++) {
                for (len_t i = 0; i < n1[ir]; i++) {
                    F2(ir, i, j) -= mg->GetP2_f(j) * nu_s_f2[ir][j*n1[ir]+i];
                }
            }
        }    
    }
}


/**
 * Sets the Jacobian matrix for the specified block
 * in the given matrix.
 * NOTE: This routine assumes that the advection coefficients
 * are independent of all other unknown quantities (solved
 * for at the same time).
 *
 * uqtyId:  ID of the unknown quantity which the term
 *          is applied to (block row).
 * derivId: ID of the quantity with respect to which the
 *          derivative is to be evaluated.
 * mat:     Jacobian matrix block to populate.
 * x:       Value of the unknown quantity.
 * 
 */
void SlowingDownTerm::SetJacobianBlock(
    const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t* x
) {

    if (uqtyId == derivId)
        this->SetMatrixElements(jac, nullptr);
        
   
   /**
    * Check if derivId is one of the id's that contributes 
    * to this advection coefficient 
    */
   if( (derivId == id_ncold) || (derivId == id_ni) ){
        
        // The number of subquantities in derivId (i.e. number of ion species)
        len_t nMultiples = 1*(derivId == id_ncold) + nzs*(derivId == id_ni);

        // Allocate
        real_t **dfr = new real_t*[(nr+1)*nMultiples]; 
        real_t **df1 = new real_t*[nr*nMultiples];
        real_t **df2 = new real_t*[nr*nMultiples];

        for(len_t ir = 0; ir<nr; ir++){
            for(len_t n = 0; n<nMultiples; n++){        
                df1[n*nr+ir] = new real_t[(n1[ir]+1)*n2[ir]];
                df2[n*nr+ir] = new real_t[n1[ir]*(n2[ir]+1)];
            }
        }
        real_t *JacobianColumn = new real_t[grid->GetNCells()];

        // XXX: assume same n1, n2 at all radii for the (all-zero) Fr term
        for(len_t ir = 0; ir<nr+1; ir++)
            for(len_t n = 0; n<nMultiples; n++) 
                dfr[n*(nr+1)+ir] = new real_t[n1[0]*n2[0]];

        // Set partial advection coefficients for this advection term 
        GetPartialAdvectionTerm(derivId,dfr, df1,df2, nMultiples);

        for(len_t n=0; n<nMultiples; n++){

            // Initialize the jacobian column to 0
            len_t offset = 0; 
            for(len_t ir=0; ir<nr; ir++){
                for (len_t j = 0; j < n2[ir]; j++) 
                    for (len_t i = 0; i < n1[ir]; i++) 
                        JacobianColumn[offset + n1[ir]*j + i] = 0;

                offset += n1[ir]*n2[ir];
            }

            // Set one subquantity of the jacobian matrix 
            SetVectorElements(JacobianColumn, x, dfr+n*nr, df1+n*nr, df2+n*nr);
            offset = 0;
            for(len_t ir=0; ir<nr; ir++){
                for (len_t j = 0; j < n2[ir]; j++) 
                    for (len_t i = 0; i < n1[ir]; i++) 
                        jac->SetElement(offset + n1[ir]*j + i, n*nr+ir, JacobianColumn[offset + n1[ir]*j + i]); 

                offset += n1[ir]*n2[ir];
            }
        }
   

        // Deallocate 
        for(len_t ir = 0; ir<nr; ir++){
            for(len_t n = 0; n<nMultiples; n++){        
                delete [] df1[n*nr+ir];
                delete [] df2[n*nr+ir];
            }
        }
        for(len_t ir = 0; ir<nr+1; ir++)
            for(len_t n = 0; n<nMultiples; n++)
                delete [] dfr[n*(nr+1)+ir];
            
        delete [] JacobianColumn;
        delete [] df1;
        delete [] df2;
        delete [] dfr;        
    }
}


// Return jacobian of the advection coefficient for this advection term
void SlowingDownTerm::GetPartialAdvectionTerm(len_t derivId, real_t **&dfr, real_t **&df1, real_t **&df2, len_t nMultiples){
    const len_t nr = grid->GetNr();
 
    // XXX: assume same n1, n2 on all radii
    for(len_t n=0; n<nMultiples; n++){
        for(len_t ir = 0; ir<nr+1; ir++){
            for (len_t j = 0; j < n2[0]; j++) {
                for (len_t i = 0; i < n1[0]; i++) {
                    dfr[n*(nr+1)+ir][j*n1[0]+i] = 0;
                }
            }
        }
    }

    // Get jacobian of the collision frequency
    const real_t *dNuS_f1 = nuS->GetUnknownPartialContribution(derivId, FVM::FLUXGRIDTYPE_P1);
    const real_t *dNuS_f2 = nuS->GetUnknownPartialContribution(derivId, FVM::FLUXGRIDTYPE_P2);
  
    bool gridtypePXI        = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI);
    bool gridtypePPARPPERP  = (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PPARPPERP);

    len_t offset1 = 0;
    len_t offset2 = 0;

    // Set partial advection coefficients
    for(len_t n=0; n<nMultiples; n++){
        for (len_t ir = 0; ir < nr; ir++) {
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = n1[ir];
        const len_t np2 = n2[ir];
            if (gridtypePXI || gridtypePPARPPERP) {
                for (len_t j = 0; j < np2; j++) {
                    for (len_t i = 0; i < np1+1; i++) {
                        df1[nr*n+ir][j*(n1[ir]+1)+i] = -mg->GetP1_f(i) * dNuS_f1[offset1 + j*(np1+1) + i];
                    }
                }
            }

            if (gridtypePPARPPERP) {
                for (len_t j = 0; j < np2+1; j++) {
                    for (len_t i = 0; i < np1; i++) {
                        df2[nr*n+ir][j*n1[ir]+i] = -mg->GetP2_f(j) * dNuS_f2[offset2 + j*np1 + i];
                    }
                }
            } else if (gridtypePXI) {
                for (len_t j = 0; j < np2+1; j++) {
                    for (len_t i = 0; i < np1; i++) {
                        df2[nr*n+ir][j*n1[ir]+i] = 0;
                    }
                }
            }
            offset1 += (np1+1)*np2;
            offset2 += np1*(np2+1);
        }
    }
}