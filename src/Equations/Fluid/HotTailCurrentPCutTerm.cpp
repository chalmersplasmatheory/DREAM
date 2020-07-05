
#include "DREAM/Equations/Fluid/HotTailCurrentPCutTerm.hpp"
#include <limits>
using namespace DREAM;


HotTailCurrentPCutTerm::HotTailCurrentPCutTerm(
            FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, 
            FVM::UnknownQuantityHandler *u, PitchScatterFrequency *nuD
) : EquationTerm(fluidGrid), fluidGrid(fluidGrid), hottailGrid(hottailGrid), unknowns(u), nuD(nuD) {}

HotTailCurrentPCutTerm::~HotTailCurrentPCutTerm(){
    Deallocate();
}


void HotTailCurrentPCutTerm::Rebuild(const real_t,const real_t, FVM::UnknownQuantityHandler* unknowns) {
    id_fhot = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
    id_jhot = unknowns->GetUnknownID(OptionConstants::UQTY_J_HOT);
    id_pcut = unknowns->GetUnknownID(OptionConstants::UQTY_J_HOT_P_CUT);
    id_Eterm =  unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    id_ni = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    const real_t *pcut = unknowns->GetUnknownData(id_pcut);
//    const real_t *Eterm = unknowns->GetUnknownData(id_Eterm);

    if(!hasBeenInitialised)
        hasBeenInitialised = GridRebuilt();

    for(len_t ir = 0; ir<nr; ir++){
        FVM::MomentumGrid *mg = hottailGrid->GetMomentumGrid(ir);
        const real_t *p_f = mg->GetP1_f();
        const real_t *p   = mg->GetP1();
        len_t np = mg->GetNp1();
        for(len_t i=0; i<np; i++){
            iCut[ir] = i;
            if( (pcut[ir] < p_f[i+1]) && (pcut[ir] > p_f[i]) )
                break;
        }        

        if(iCut[ir]==0)
            dp[ir] = 2*p[0];
        else if(iCut[ir]==np-1)
            dp[ir] = p_f[np]-p[np-1];
        else
            dp[ir] = p[iCut[ir]+1]- p[iCut[ir]-1];

        GOverH[ir] = EvalGOverH(ir,-1);
    }
}

/**
 * Method that is called whenever the grid is rebuilt. 
 * Allocates memory. 
 */
bool HotTailCurrentPCutTerm::GridRebuilt() {
    Deallocate();
    nr = fluidGrid->GetNr();
    dp = new real_t[nr];
    iCut = new len_t[nr];
    GOverH = new real_t[nr];
//    dGOverH = new real_t[nr];
    return true;
}

void HotTailCurrentPCutTerm::Deallocate(){
    if(dp == nullptr)
        return;

    delete [] dp;
    delete [] iCut;
    delete [] GOverH;
//    delete [] dGOverH;
    
}


real_t HotTailCurrentPCutTerm::EvalGOverH(len_t ir, len_t derivId, len_t n){
    real_t pcut = unknowns->GetUnknownData(id_pcut)[ir];
    const real_t EffPass = fluidGrid->GetRadialGrid()->GetEffPassFrac(ir);
    const real_t Bavg = fluidGrid->GetRadialGrid()->GetFSA_B2(ir);
    real_t constTerm = Constants::ec /(Constants::me * Constants::c)
            * EffPass / (3*sqrt(Bavg) );
    real_t nud = nuD->evaluateAtP(ir,pcut);
    real_t oneOverNuD = 1/nud;
    real_t E = unknowns->GetUnknownData(id_Eterm)[ir];
    if(derivId==id_Eterm){
        E = 1;
    } else if((derivId==id_ncold) || (derivId==id_ni)) {
        oneOverNuD = -nuD->evaluatePartialAtP(ir,pcut,derivId,n) /(nud*nud);
    } else if(derivId==id_pcut){
        real_t eps = std::numeric_limits<real_t>::epsilon();
        real_t dp = pcut*sqrt(eps);
        oneOverNuD = (1/nuD->evaluateAtP(ir,pcut+dp) - 1/nud)/dp;
    }
    return constTerm * E * oneOverNuD;
}

/**
 * Set the jacobian elements for this term. 
 *
 * derivId: Unknown ID of derivative with respect to which differentiation
 *          should be done.
 * unknId:  ID of the unknown to differentiate.
 * jac:     Jacobian matrix to set elements of.
 * x:       Value of the unknown quantity.
 */
void HotTailCurrentPCutTerm::SetJacobianBlock(
    const len_t unknId, const len_t derivId, FVM::Matrix *jac, const real_t* f
) {
    if (derivId == unknId)
        this->SetMatrixElements(jac,nullptr);

    if( !((derivId==id_pcut) || (derivId==id_Eterm) || (derivId==id_ni) || (derivId==id_ncold)) )
        return;
    
    len_t offset_n = 0;
    for(len_t n=0; n<unknowns->GetUnknown(derivId)->NumberOfMultiples(); n++){    
        len_t offset = 0;
        for(len_t ir=0;ir<nr;ir++){        
            len_t np = hottailGrid->GetMomentumGrid(ir)->GetNp1();
            len_t ind = offset + iCut[ir];
            real_t dGOverH = EvalGOverH(ir,derivId,n);
            real_t fPlus, fMinus;
            if(iCut[ir]==0)
                fMinus = f[ind];
            else
                fMinus = f[ind-1];
            if(iCut[ir]==np-1)
                fPlus = 0;
            else
                fPlus = f[ind+1];
            
            real_t F = dGOverH * (fPlus - fMinus)/dp[ir];
            jac->SetElement(ir,n*nr + ir, F);

            offset += np;
        }
        
        offset_n += hottailGrid->GetNCells();
    }
}


/**
 * Set the elements of the linear operator matrix corresponding to
 * this operator.
 *
 * mat: Linear operator matrix to set elements of.
 * rhs: Equation right-hand-side.
 */
void HotTailCurrentPCutTerm::SetMatrixElements(FVM::Matrix *mat, real_t*){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        len_t np = hottailGrid->GetMomentumGrid(ir)->GetNp1();
        len_t ind = offset + iCut[ir];
        len_t indMinus = ind-1;
        if(iCut[ir] < np-1) // no contribution on upper boundary
            mat->SetElement(ir, ind+1, GOverH[ir]/dp[ir]);
        if(iCut[ir]==0) // neumann condition at i=0
            indMinus = 0;
        mat->SetElement(ir, indMinus, -GOverH[ir]/dp[ir]);
        mat->SetElement(ir, ind, 1);
        offset += np;
    }
    
}


void HotTailCurrentPCutTerm::SetVectorElements(real_t *vec, const real_t *f){

    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        len_t np = hottailGrid->GetMomentumGrid(ir)->GetNp1();
        len_t ind = offset + iCut[ir];
        real_t fPlus, fMinus;
        if(iCut[ir]==0)
            fMinus = f[ind];
        else
            fMinus = f[ind-1];
        if(iCut[ir]==np-1)
            fPlus = 0;
        else
            fPlus = f[ind+1];
        
        vec[ir] += GOverH[ir] * (fPlus - fMinus)/dp[ir]  + f[ind];    
        offset += np;
    }

}
