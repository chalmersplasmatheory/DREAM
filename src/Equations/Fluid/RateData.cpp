#include "DREAM/Equations/Fluid/RateData.hpp"
//would be nice if this is also generated from the tools

  namespace DREAM {

    const MolecularRatePairDefinition molecularRatePairDefinitions[] = {
      /*The formula is as follows: 
      Reactant1 + Reactant2 -> 
      Product1 + Product2
      Charge exchange, 
      Dissociation,
      Dissociative recombination rates 
      */

      {
          {"Ar", 0}, {"Ar", 1}, 
          {"Ar", 1}, {"Ar", 0},
          {"Ar_Ar_charge_exchange"},      
          {nullptr}, 
          {nullptr}  
      },
      {
          {"Ar", 0}, {"Ar2", 1}, {"Ar2", 1}, {"Ar", 0},
          {"Ar_Ar2_charge_exchange"},  
          {nullptr}, 
          {nullptr}      
      },
      {
          {"D", 0}, {"D", 1}, {"D", 1}, {"D", 0},
          {"D_D_charge_exchange"},   
          {nullptr}, 
          {nullptr}   
      },
      {
          {"D2", 0}, {"D2", 1}, {"D2", 1}, {"D2", 0},
          {"D2_D2_charge_exchange"},    
          {nullptr},
          {nullptr}
      },
      {
          {"D2", 1}, {"D", 0},{"D2", 0}, {"D", 1},
          {"D2_D_charge_exchange"}, 
          {nullptr},
          {nullptr}       
      },

      {
          {"Ar", 1}, {"D", 0},{"Ar", 0}, {"D", 1},
          {"Ar_D_charge_exchange"},    
          {nullptr}, 
          {nullptr}      
      },

      {
          {"Ar", 2}, {"D", 0},{"Ar", 1}, {"D", 1},
          {"Ar2_D_charge_exchange"},     
          {nullptr},
          {nullptr}      
      },     
      {
          {"Ar", 3}, {"D", 0},{"Ar", 2}, {"D", 1},
          {"Ar3_D_charge_exchange"}, 
          {nullptr}, 
          {nullptr}      
      },

      {
          {"Ar", 4}, {"D", 0},{"Ar", 3}, {"D", 1},
          {"Ar4_D_charge_exchange"}, 
          {nullptr}, 
          {nullptr}      
      },

      {
          {"e", 0}, {"D2", 0},{"D", 0}, {"D", 0},
          {nullptr}, 
          {"e_D2_dissociation"},
          {nullptr}      
      },

      {
          {"e", 0}, {"D2", 0},{"D", 0}, {"D", 1},
          {nullptr}, 
          {nullptr},
          {"e_D2_dissociative_ionization"}      
      },
      {
          {"D", 1}, {"D2", 0},{"D", 0}, {"D2", 1}, 
          {"D_D2_charge_exchange"} , 
          {nullptr} ,
          {nullptr}      
      },
      {
          {"D2", 1}, {"D2", 0},{"D3", 1}, {"D", 0}, 
          {nullptr}, 
          {"D2_D2_D3_dissociation"}, 
          {nullptr}      
      },
      {
          {"Ar", 1}, {"D2", 0},{"Ar", 0}, {"D2", 1}, 
          {"Ar_D2_charge_exchange"},
          {nullptr}, 
          {nullptr}      
      },

      {
          {"Ar", 1}, {"D2", 0},{"ArD", 1}, {"D", 0}, 
          {nullptr},
          {"Ar_D2_dissociation"}, 
          {nullptr}      
      },
      {
          {"Ar", 2}, {"D2", 0},{"Ar", 1}, {"D2", 1}, 
          {"Ar2_D2_charge_exchange"},
          {nullptr}, 
          {nullptr}      
      },

      {
          {"Ar", 3}, {"D2", 0},{"Ar", 2}, {"D2", 1}, 
          {"Ar3_D2_charge_exchange"},
          {nullptr}, 
          {nullptr}      
      },
      {
          {"Ar", 4}, {"D2", 0},{"Ar", 3}, {"D2", 1}, 
          {"Ar4_D2_charge_exchange"},
          {nullptr}, 
          {nullptr}      
      },

      {
          {"ArD", 1}, {"D2", 0},{"Ar", 0}, {"D3", 1}, 
          {nullptr}, //TODO
          {"ArD_D2_dissociation"}, 
          {nullptr}      
      },

      {
          {"D2", 1}, {"Ar", 0},{"D", 0}, {"ArD", 1}, 
          {nullptr}, 
          {"D2_Ar_dissociation"}, 
          {nullptr}      
      },
      {
          {"D3", 1}, {"Ar", 0},{"D2", 0}, {"ArD", 1}, 
          {nullptr}, 
          {"D3_Ar_dissociation"}, 
          {nullptr}      
      },

  };

  const len_t molecularRatePairDefinitionCount =
      sizeof(molecularRatePairDefinitions) /
      sizeof(molecularRatePairDefinitions[0]);
  


}