#ifndef _DREAM_EQUATIONS_MOLECULE_HANDLER_HPP
#define _DREAM_EQUATIONS_MOLECULE_HANDLER_HPP

#include <vector>
#include <string>
#include "FVM/config.h"


namespace DREAM {

class MoleculeHandler {
public:
    bool IsMolecule(const std::string&) const;
    real_t GetMass(const std::string&) const;

    MoleculeHandler() = default;
    ~MoleculeHandler() = default;

  };



} // namespace DREAM

#endif // _DREAM_EQUATIONS_MOLECULE_HANDLER_HPP