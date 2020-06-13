
#include "DREAM/Equations/Fluid/BindingEnergyTerm.hpp"
#include "DREAM/NIST.hpp"


using namespace DREAM;


// number of species in dataset
/*const len_t BindingEnergyTerm::nZ = 2; 
// total number of data points = sum(dataZ)
const len_t BindingEnergyTerm::nData = 19; 
// species in dataset
const real_t BindingEnergyTerm::dataZ[nZ] = {1,18}; */
/**
 * Ordered list of total binding energies for ions in eV, 
 * starting from Z0=0 to Z0=Z-1 (where Z0=Z is always 0, 
 * handled by GetBindingEnergy.)
 * Data retrieved from NIST website ("Total Binding Energy"):
 * https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
 */
//const real_t BindingEnergyTerm::dataIbind[nData] = {/*H*/ 13.59843 /*end H*/,/*Ar*/ 14400.8, 14385.0,14357.4, 14316.7, 14257.1, 14182.2, 14090.9, 13966.6, 13823.1, 13400.5, 12920.8, 12380.3, 11761.3, 11075.9, 10320.7, 9465.264, 8546.8885, 4426.2228 /*end Ar*/};


/**
 * Constructor.
 */
BindingEnergyTerm::BindingEnergyTerm(FVM::Grid* g, IonHandler *ionHandler, NIST *nist)
    : FVM::DiagonalLinearTerm(g), ionHandler(ionHandler), nist(nist) {
    /*if(!DataSetIsValid())
        throw FVM::FVMException("Binding energy dataset is invalid.");*/
}

/** 
 * Sets weights to the binding potential, which is _minus_
 * the binding energies since they are defined as positive.
 */
void BindingEnergyTerm::SetWeights(){
    len_t N = grid->GetNCells();
    len_t nZ = ionHandler->GetNZ();
    const len_t *Zs = ionHandler->GetZs();

    for(len_t iz = 0; iz<nZ; iz++){
        for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
            len_t n  = ionHandler->GetIndex(iz,Z0);
            real_t w = nist->GetBindingEnergy(Zs[iz],Z0);

            for(len_t i = 0; i<N; i++)
                weights[n*N+i] = -w;
        }
    }
}

/*real_t BindingEnergyTerm::GetBindingEnergy(len_t Zin, len_t Z0in){
    if(Z0in==Zin)
        return 0;
    len_t offset = 0;
    for(len_t iZ=0; iZ<nZ; iZ++){
        if(dataZ[iZ]!=Zin){ 
            // if not the requested species, skip over to next Z
            offset += dataZ[iZ];
            continue;
        } 
        for(len_t Z0 = 0; Z0<dataZ[iZ]; Z0++){
            if(Z0==Z0in)
                return dataIbind[offset];
            offset++;
        }
    }
    throw FVM::FVMException("Binding energy data is missing for Z=%u and Z0=%u in BindingEnergyTerm",Zin, Z0in);
    return 0; 
}*/


/**
 * Verifies that the format of the data is correct.
 */
/*bool BindingEnergyTerm::DataSetIsValid(){
    bool success = true; 

    // Verify that data dimensions match
    len_t nData_t = 0;
    for(len_t iz=0; iz<nZ; iz++){
        nData_t +=dataZ[iz];
    }
    if(nData_t != nData)
        success = false;

    // Verify that total binding energies decrease with increasing Z0 
    len_t offset = 0;
    for(len_t iZ=0; iZ<nZ; iZ++){
        for(len_t Z0 = 0; Z0<(dataZ[iZ]-1); Z0++){
            if(dataIbind[offset+1] > dataIbind[offset])
                success = false;
            offset++;
        }
        offset++;
    }
    

    return success;
}*/

/**
 * Set matrix for all nnz ion species
 */
void BindingEnergyTerm::SetMatrixElements(FVM::Matrix* mat, real_t*){  
    len_t N   = grid->GetNCells();
    len_t nnz = ionHandler->GetNzs();

    for(len_t n=0; n<nnz;n++)
        for(len_t i=0; i<N; i++)
            mat->SetElement(i,n*N+i,weights[n*N+i]);
}

/**
 * Set vector for all nnz ion species
 */
void BindingEnergyTerm::SetVectorElements(real_t* vec, const real_t* x){ 
    len_t N   = grid->GetNCells();
    len_t nnz = ionHandler->GetNzs();

    for(len_t i=0; i<N*nnz; i++)
        vec[i] += weights[i] * x[i];
}

