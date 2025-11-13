#ifndef _DREAM_QUANTITY_DATA_HPP
#define _DREAM_QUANTITY_DATA_HPP

#include <petsc.h>
#include <vector>
#include <softlib/SFile.h>
#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Grid/fluxGridType.enum.hpp"

namespace DREAM::FVM {
    class QuantityData {
    private:
        Grid *grid;
        std::vector<real_t> times;
        std::vector<real_t*> store;
        enum FVM::fluxGridType fluxGridType = FLUXGRIDTYPE_DISTRIBUTION;

        len_t nMultiples=1;
        len_t nElements=0;
        // Data in current step
        real_t *data=nullptr;
        // Data in previous time step(s) (even if step was not saved to 'store')
        // (DREAM always needs access to the previous and current time steps in
        // order to calculate time derivatives. In order to be able to roll back
        // solutions to a previous state, we can increase the number of old solutions
        // stored here. This is necessary for the adaptive time stepper, but also
        // implies a higher memory consumption)
        real_t **olddata=nullptr;
        real_t *oldtime = nullptr;
        const len_t N_SAVE_OLD_STEPS = 4;   // Can roll back N-1 steps (TimeStepperAdaptive needs N >= 3 (so that we can also restore the "initial" time derivative))
        len_t nOldSaved = 0;      // Number of old steps currently stored

        // Data from time step before the previous (even in step was not saved to 'store')
        // (this variable is one time step older than 'olddata' and is used when
        // rolling back saved steps)

        // Vector used for addressing PETSc vectors
        PetscInt *idxVec = nullptr;

        bool hasChanged = true;

        void AllocateData();

        void SaveSFile_internal(SFile*, const std::string& name, const std::string&, const std::string&, bool saveMeta, std::vector<real_t>&, std::vector<real_t*>&);

    public:
        QuantityData(
            FVM::Grid*, const len_t nMultiples=1,
            enum FVM::fluxGridType fgt=FLUXGRIDTYPE_DISTRIBUTION
        );
        ~QuantityData();

        real_t *Get() { return this->data; }
        //real_t *GetPrevious() { return this->store.back(); }
        real_t *GetPrevious() { return this->olddata[0]; }
        real_t GetPreviousTime() { return this->oldtime[0]; }
        real_t *GetInitialData() { return this->store.front(); }
        len_t Size() { return this->nElements; }
		len_t GetNumberOfSavedSteps() { return this->store.size(); }

        /**
         * Returns 'true' if the data stored by this quantity
         * was changed in the last call to 'Store()'.
         */
        bool HasChanged() const { return this->hasChanged; }
        bool HasInitialValue() const { return (this->store.size()>=1); }

        len_t GetNOldSaved() const { return this->nOldSaved; }

        bool CanRollbackSaveStep() const;
        void RollbackSaveStep();
        void SaveStep(const real_t, bool);
        void Store(Vec&, const len_t, bool mayBeConstant=false);
        void Store(const real_t*, const len_t offset=0, bool mayBeConstant=false);
        void Store(const int_t*, const len_t offset=0, bool mayBeConstant=false);
        void Store(const len_t, const len_t, const real_t *const*, bool mayBeConstant=false);
        real_t *StoreEmpty();
		void RestoreValue();

        void SaveSFile(SFile*, const std::string& name, const std::string& path="", const std::string& desc="", bool saveMeta=false);
		void SaveSFileCurrent(SFile*, const std::string& name, const std::string& path="", const std::string& desc="", bool saveMeta=false);

        void SetInitialValue(const real_t*, const real_t t0=0);
    };
}

#endif/*_DREAM_QUANTITY_DATA_HPP*/
