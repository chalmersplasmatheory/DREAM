#ifndef _DREAM_FVM_UNKNOWN_QUANTITY_HANDLER_HPP
#define _DREAM_FVM_UNKNOWN_QUANTITY_HANDLER_HPP

namespace DREAM::FVM { class UnknownQuantityHandler; }

#include <string>
#include <vector>
#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantity.hpp"

namespace DREAM::FVM {
    class UnknownQuantityHandler {
    private:
        std::vector<UnknownQuantity*> unknowns;

    public:
        UnknownQuantityHandler();
        ~UnknownQuantityHandler();

        UnknownQuantity *GetUnknown(const len_t i) { return unknowns.at(i); }
        len_t GetUnknownID(const std::string&);
        len_t GetNUnknowns() const { return this->unknowns.size(); }
        len_t Size() const { return GetNUnknowns(); }

        const real_t *GetLongVector(const std::vector<len_t>& nontrivials, real_t *vec=nullptr);
        const real_t *GetLongVector(const len_t, const len_t*, real_t *vec=nullptr);
        const real_t *GetLongVectorPrevious(const std::vector<len_t>& nontrivials, real_t *vec=nullptr);
        const real_t *GetLongVectorPrevious(const len_t, const len_t*, real_t *vec=nullptr);
        const len_t GetLongVectorSize(const std::vector<len_t>& nontrivials);
        const len_t GetLongVectorSize(const len_t, const len_t*);

        const real_t *GetLongVectorAll(real_t *vec=nullptr);
        const len_t GetLongVectorSizeAll();

        real_t GetUnknownDataPreviousTime(const len_t);
        real_t *GetUnknownData(const len_t);
        real_t *GetUnknownData(const std::string&);
        real_t *GetUnknownDataPrevious(const len_t);
        real_t *GetUnknownDataPrevious(const std::string&);
        real_t *GetUnknownInitialData(const len_t);
        const std::vector<UnknownQuantity*>& GetUnknowns() const { return this->unknowns; }

        bool HasChanged(const len_t id) const { return unknowns[id]->HasChanged(); }
        bool HasInitialValue(const len_t id) const { return unknowns[id]->HasInitialValue(); }
        bool HasUnknown(const std::string&);
        len_t InsertUnknown(const std::string&, const std::string&, Grid*, const len_t nMultiples=1);

        void Store(std::vector<len_t>&, Vec&, bool mayBeConstant=false);
        void Store(std::vector<len_t>&, const real_t*, bool mayBeConstant=false);
        void Store(const len_t id, Vec& v, const len_t offs, bool mayBeConstant=false) { unknowns[id]->Store(v, offs, mayBeConstant); }
        void Store(const len_t id, const real_t *v, const len_t offs=0, bool mayBeConstant=false) { unknowns[id]->Store(v, offs, mayBeConstant); }
        void SetFromLongVectorAll(const real_t*, bool mayBeConstant=false);
		void RestoreSolution(std::vector<len_t>&);

        void RollbackSaveStep();
        void SaveStep(const real_t t, bool trueSave);

        void SaveSFile(const std::string& filename, bool saveMeta=false);
        void SaveSFile(SFile*, const std::string& path="", bool saveMeta=false);
        void SaveSFileCurrent(const std::string& filename, bool saveMeta=false);
        void SaveSFileCurrent(SFile*, const std::string& path="", bool saveMeta=false);

        void SetInitialValue(const std::string&, const real_t*, const real_t t0=0);
        void SetInitialValue(const len_t, const real_t*, const real_t t0=0);

        UnknownQuantity *operator[](const len_t i) { return GetUnknown(i); }
    };
}

#endif/*_DREAM_FVM_UNKNOWN_QUANTITY_HANDLER_HPP*/
