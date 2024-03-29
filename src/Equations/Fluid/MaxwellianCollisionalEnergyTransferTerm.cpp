#include "DREAM/Equations/Fluid/MaxwellianCollisionalEnergyTransferTerm.hpp"
#include "DREAM/Constants.hpp"

/**
 * Implementation of an equation terms which represents the rate of energy transfer (in SI)
 * between two maxwellians of densities n_i, n_j and heat content W_i, W_j with W = 1.5nT
 * representing particles with masses m_i, m_j and charge numbers Z_i and Z_j.
 * 
 * It can describe both ion-ion and ion-electron heat transfer, only differing in 
 * definition by the Coulomb logarithm used, controlled with the `bool isEI` flag.
 * 
 * One instance of this class represents the interaction between just two species.
 */

using namespace DREAM;

/**
 * Constructor.
 * The unknown id's defines which quantities are targets and which are "incident" 
 * (should be either ions or electrons). The offsets should be 0 for electrons and
 * the ion index for ions. 
 */
MaxwellianCollisionalEnergyTransferTerm::MaxwellianCollisionalEnergyTransferTerm(
            FVM::Grid *g, 
            const len_t index_i, bool isIon_i,
            const len_t index_j, bool isIon_j, 
            FVM::UnknownQuantityHandler *u, CoulombLogarithm *lnL, IonHandler *ionHandler, real_t scaleFactor
) : FVM::EquationTerm(g), index_i(index_i), index_j(index_j), isIon_i(isIon_i), isIon_j(isIon_j),
    unknowns(u), lnLambda(lnL), ionHandler(ionHandler)
{
    SetName("MaxwellianCollisionalEnergyTransferTerm");

    if(isIon_i){
        this->mi = ionHandler->GetIonSpeciesMass(index_i);
        this->Zi = ionHandler->GetZ(index_i);
    } else {
        this->mi = Constants::me;
        this->Zi = 1;
    } if(isIon_j){
        this->mj = ionHandler->GetIonSpeciesMass(index_j);
        this->Zj = ionHandler->GetZ(index_j);
    } else { 
        this->mj = Constants::me;
        this->Zj = 1;
    }
    real_t ec = Constants::ec;
    real_t eps0 = Constants::eps0;
    real_t ec2eps = ec*ec / eps0;
    this->constPreFactor = scaleFactor*sqrt(3.0/M_PI)*ec2eps*ec2eps*sqrt(mi*mj) / (4*M_PI);
    this->isEI = !(isIon_i && isIon_j);
    this->lnLambda_settings = new CollisionQuantity::collqty_settings;
    this->lnLambda_settings->lnL_type = isEI ? 
        OptionConstants::COLLQTY_LNLAMBDA_THERMAL : OptionConstants::COLLQTY_LNLAMBDA_ION_ION;

    this->id_ncold = u->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_Wcold = u->GetUnknownID(OptionConstants::UQTY_W_COLD);
    this->id_Ni    = u->GetUnknownID(OptionConstants::UQTY_NI_DENS);
    this->id_Wi    = u->GetUnknownID(OptionConstants::UQTY_WI_ENER);
    this->id_ions  = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->id_Tcold = u->GetUnknownID(OptionConstants::UQTY_T_COLD);

    AddUnknownForJacobian(u, id_ncold);
    AddUnknownForJacobian(u, id_Ni);
    AddUnknownForJacobian(u, id_Wi);
    AddUnknownForJacobian(u, id_Wcold);
    AddUnknownForJacobian(u, id_ions);
    AddUnknownForJacobian(u, id_Tcold);    
}


/**
 * Set vector elements of this equation term. The implementation mirrors the analytic
 * form that is given in doc/notes/theory, where we have rewritten it in terms of the
 * total ion densities Ni and heat Wi.
 */
void MaxwellianCollisionalEnergyTransferTerm::SetVectorElements(real_t *vec, const real_t*){
    const real_t *lnL = isEI ? lnLambda->GetLnLambdaT() : lnLambda->GetLnLambdaII();
    for(len_t ir=0; ir<nr; ir++){
        real_t nZ2_i,nZ2_j,ni, nj, Wi, Wj;
        GetParametersForSpecies(ir, index_i, isIon_i, ni, Wi, nZ2_i);
        GetParametersForSpecies(ir, index_j, isIon_j, nj, Wj, nZ2_j);

        real_t njWi = nj*Wi;
        real_t niWj = ni*Wj;
        if(ni==0 || nj==0 || niWj==njWi)
            continue;
        real_t pre = sqrt(ni*nj)*nZ2_i*nZ2_j;
        real_t up = niWj - njWi;
        real_t down = mj*njWi + mi*niWj;
        vec[index_i*nr + ir] += lnL[ir] * constPreFactor * pre * up / (down*sqrt(down));
    }
}


/**
 * Sets the jacobian block of this equation term
 */
bool MaxwellianCollisionalEnergyTransferTerm::SetJacobianBlock(const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t*){
    if( !HasJacobianContribution(derivId) ) 
        return false;

    const real_t *lnL = isEI ? lnLambda->GetLnLambdaT() : lnLambda->GetLnLambdaII();
    
    for(len_t ir=0; ir<nr; ir++){
        real_t nZ2_i,nZ2_j,ni, nj, Wi, Wj;
        GetParametersForSpecies(ir, index_i, isIon_i, ni, Wi, nZ2_i);
        GetParametersForSpecies(ir, index_j, isIon_j, nj, Wj, nZ2_j);

        real_t njWi = nj*Wi;
        real_t niWj = ni*Wj;
        if(ni==0 || nj==0 || niWj==njWi) // last one (Ti=Tj) deals with Wi=Wj=0 which is singular below
            continue;
        real_t pre = sqrt(ni*nj)*nZ2_i*nZ2_j;
        real_t up = niWj - njWi;
        real_t down = mj*njWi + mi*niWj;
        real_t vec = lnL[ir] * constPreFactor * pre * up / (down*sqrt(down));

        len_t ii = index_i*nr + ir;
        len_t jj = index_j*nr + ir;
        real_t val1 = 0, val2 = 0;
        if( (isIon_i &&  derivId==id_Ni) || (!isIon_i && derivId==id_ncold) ) 
            val1 += vec*( 0.5/ni + Wj*( 1.0/up - 1.5*mi/down) );
        if( (isIon_j &&  derivId==id_Ni) || (!isIon_j && derivId==id_ncold) ) 
            val2 += vec*( 0.5/nj + Wi*(-1.0/up - 1.5*mj/down) );
        if( (isIon_i&&derivId==id_Wi) || (!isIon_i&&derivId==id_Wcold) )
            val1 += vec*nj*(-1.0/up - 1.5*mj/down );
        if( (isIon_j&&derivId==id_Wi) || (!isIon_j&&derivId==id_Wcold) )
            val2 += vec*ni*( 1.0/up - 1.5*mi/down );
        jac->SetElement(ii, ii, val1);
        jac->SetElement(ii, jj, val2);

        // handle the nZ2 term:
        if(isIon_i && derivId==id_ions){
            real_t preOverNZ2 = sqrt(ni*nj)*nZ2_j; // rebuild 'vec' without the nZ2_i factor
            real_t vecOverNZ2 = lnL[ir] * constPreFactor * preOverNZ2 * up / (down*sqrt(down));
            for(len_t Z0=0; Z0<=Zi; Z0++){
                len_t indZ = ionHandler->GetIndex(index_i, Z0);
                jac->SetElement(ii, indZ*nr + ir, vecOverNZ2 * Z0*Z0);
            }
        } else if(!isIon_i && derivId == id_ncold)
            jac->SetElement(ii,ir, vec/nZ2_i);
        if(isIon_j && derivId==id_ions){
            real_t preOverNZ2 = sqrt(ni*nj)*nZ2_i; // rebuild 'vec' without the nZ2_j factor
            real_t vecOverNZ2 = lnL[ir] * constPreFactor * preOverNZ2 * up / (down*sqrt(down));
            for(len_t Z0=0; Z0<=Zj; Z0++){
                len_t indZ = ionHandler->GetIndex(index_j, Z0);
                jac->SetElement(ii, indZ*nr + ir, vecOverNZ2 * Z0*Z0);
            }
        } else if(!isIon_j && derivId == id_ncold)
            jac->SetElement(ii,ir, vec/nZ2_j);
        
        // Below: lnLambda derivatives
        if(derivId==id_Tcold || derivId==id_ions){
            len_t nMultiple = unknowns->GetUnknown(derivId)->NumberOfMultiples();
            for(len_t n=0; n<nMultiple; n++)
                jac->SetElement(ii,n*nr+ir, vec/lnL[ir]*lnLambda->evaluatePartialAtP(ir,0,derivId,n,lnLambda_settings) );
        }
    }

    return true;
}

/**
 * Evaluates the density and heat content at radial grid point 'ir' for the specified species 
 * (ions if isIon and otherwise electrons). For ionz, nZ2 is defined as the weighted sum 
 *  nZ2_i = sum_j n_i^(j)*Z_0j^2
 * which is essentially 'Zeff' of species 'i' (multiplied by N_i). For electrons it is simply
 * taken as n_cold.
 */
void MaxwellianCollisionalEnergyTransferTerm::GetParametersForSpecies(len_t ir, len_t index,bool isIon, real_t &n, real_t &W, real_t &nZ2){
    if(isIon){
        n = unknowns->GetUnknownData(id_Ni)[nr*index + ir];
        W = unknowns->GetUnknownData(id_Wi)[nr*index + ir];
        nZ2=ionHandler->GetNZ0Z0(index, ir);
    } else {
        n = unknowns->GetUnknownData(id_ncold)[ir];
        nZ2 = n;
        W = unknowns->GetUnknownData(id_Wcold)[ir];
    }
}


