/**
 * Implementation of the object which helps map molecule densities to molecule names and charges. 
 */

#include <vector>
#include <string>
#include "DREAM/MoleculeHandler.hpp"
#include "DREAM/Constants.hpp"
#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;


struct molecule_data {
      const char *name;
      real_t mass;
  };
// in kg
static const molecule_data moleculeTable[] = {
      {"H2", 2 * Constants::mH},
      {"D2", 2 * Constants::mD},
      {"T2", 2 * Constants::mT},
      {"D3", 3 * Constants::mD},
      {"ArD", 6.97e-26},
      {"HeD",  9.99e-27},
  
  };

static const len_t nMolecules =
      sizeof(moleculeTable) / sizeof(moleculeTable[0]);

bool MoleculeHandler::IsMolecule(const std::string& name) const {
      for (len_t i = 0; i < nMolecules; i++)
          if (name == moleculeTable[i].name){
            printf("MoleculeHandler::IsMolecule: Found molecule named '%s' with mass %e.\n", name.c_str(), moleculeTable[i].mass);
              return true;
          }
      return false;
  }


real_t MoleculeHandler::GetMass(const std::string& name) const {
    for (len_t i = 0; i < nMolecules; i++)
          if (name == moleculeTable[i].name)
              return moleculeTable[i].mass;

    throw FVM::FVMException(
          "MoleculeHandler: No molecule named '%s' has been defined.",
          name.c_str()
      );
  }


