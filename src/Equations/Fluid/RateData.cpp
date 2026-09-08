#include "DREAM/Equations/Fluid/RateData.hpp"

namespace DREAM {

static const MolecularReactionSpecies Ar_Arp_charge_exchange_reactants[] = {
    {"Ar", 0, 1},
    {"Ar", 1, 1},
};

static const MolecularReactionSpecies Ar_Arp_charge_exchange_products[] = {
    {"Ar", 1, 1},
    {"Ar", 0, 1},
};

static const MolecularReactionSpecies Ar_Ar2p_charge_exchange_reactants[] = {
    {"Ar", 0, 1},
    {"Ar2", 1, 1},
};

static const MolecularReactionSpecies Ar_Ar2p_charge_exchange_products[] = {
    {"Ar", 1, 1},
    {"Ar2", 0, 1},
};

static const MolecularReactionSpecies D_Dp_charge_exchange_reactants[] = {
    {"D", 0, 1},
    {"D", 1, 1},
};

static const MolecularReactionSpecies D_Dp_charge_exchange_products[] = {
    {"D", 1, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies D_2_D_2p_charge_exchange_reactants[] = {
    {"D2", 0, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies D_2_D_2p_charge_exchange_products[] = {
    {"D2", 1, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies D_2p_D_charge_exchange_reactants[] = {
    {"D2", 1, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies D_2p_D_charge_exchange_products[] = {
    {"D2", 0, 1},
    {"D", 1, 1},
};

static const MolecularReactionSpecies Arp_D_charge_exchange_reactants[] = {
    {"Ar", 1, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies Arp_D_charge_exchange_products[] = {
    {"Ar", 0, 1},
    {"D", 1, 1},
};

static const MolecularReactionSpecies Ar2p_D_charge_exchange_reactants[] = {
    {"Ar", 2, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies Ar2p_D_charge_exchange_products[] = {
    {"Ar", 1, 1},
    {"D", 1, 1},
};

static const MolecularReactionSpecies Ar3p_D_charge_exchange_reactants[] = {
    {"Ar", 3, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies Ar3p_D_charge_exchange_products[] = {
    {"Ar", 2, 1},
    {"D", 1, 1},
};

static const MolecularReactionSpecies Ar4p_D_charge_exchange_reactants[] = {
    {"Ar", 4, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies Ar4p_D_charge_exchange_products[] = {
    {"Ar", 3, 1},
    {"D", 1, 1},
};

static const MolecularReactionSpecies D_2_dissociation_reactants[] = {
    {"e", -1, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies D_2_dissociation_products[] = {
    {"e", -1, 1},
    {"D", 1, 2},
};

static const MolecularReactionSpecies D_2_ionization_reactants[] = {
    {"e", -1, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies D_2_ionization_products[] = {
    {"e", -1, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies D_2_dissociation_2_reactants[] = {
    {"e", -1, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies D_2_dissociation_2_products[] = {
    {"D", 0, 1},
    {"D", 1, 1},
};

static const MolecularReactionSpecies Dp_D_2_charge_exchange_reactants[] = {
    {"D", 1, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies Dp_D_2_charge_exchange_products[] = {
    {"D", 0, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies D_2p_D_2_dissociation_reactants[] = {
    {"D2", 1, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies D_2p_D_2_dissociation_products[] = {
    {"D3", 1, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies Arp_D_2_charge_exchange_reactants[] = {
    {"Ar", 1, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies Arp_D_2_charge_exchange_products[] = {
    {"Ar", 0, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies Arp_D_2_dissociation_reactants[] = {
    {"Ar", 1, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies Arp_D_2_dissociation_products[] = {
    {"ArD", 1, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies Ar2p_D_2_charge_exchange_reactants[] = {
    {"Ar", 2, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies Ar2p_D_2_charge_exchange_products[] = {
    {"Ar", 1, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies Ar3p_D_2_charge_exchange_reactants[] = {
    {"Ar", 3, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies Ar3p_D_2_charge_exchange_products[] = {
    {"Ar", 2, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies Ar4p_D_2_charge_exchange_reactants[] = {
    {"Ar", 4, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies Ar4p_D_2_charge_exchange_products[] = {
    {"Ar", 3, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies ArDp_D_2_dissociation_reactants[] = {
    {"ArD", 1, 1},
    {"D2", 0, 1},
};

static const MolecularReactionSpecies ArDp_D_2_dissociation_products[] = {
    {"Ar", 0, 1},
    {"D3", 1, 1},
};

static const MolecularReactionSpecies D_2p_Ar_dissociation_reactants[] = {
    {"D2", 1, 1},
    {"Ar", 0, 1},
};

static const MolecularReactionSpecies D_2p_Ar_dissociation_products[] = {
    {"D", 0, 1},
    {"ArD", 1, 1},
};

static const MolecularReactionSpecies D_3p_Ar_dissociation_reactants[] = {
    {"D3", 1, 1},
    {"Ar", 0, 1},
};

static const MolecularReactionSpecies D_3p_Ar_dissociation_products[] = {
    {"D2", 0, 1},
    {"ArD", 1, 1},
};

static const MolecularReactionSpecies D_2p_recombination_reactants[] = {
    {"e", -1, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies D_2p_recombination_products[] = {
    {"D", 0, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies D_2p_recombination_2_reactants[] = {
    {"e", -1, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies D_2p_recombination_2_products[] = {
    {"D", 0, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies D_3p_dissociation_reactants[] = {
    {"e", -1, 1},
    {"D3", 1, 1},
};

static const MolecularReactionSpecies D_3p_dissociation_products[] = {
    {"D2", 0, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies D_3p_dissociation_2_reactants[] = {
    {"e", -1, 1},
    {"D3", 1, 1},
};

static const MolecularReactionSpecies D_3p_dissociation_2_products[] = {
    {"D", 0, 3},
};

static const MolecularReactionSpecies D_2pAr2p_dissociation_reactants[] = {
    {"D2", 1, 1},
    {"Ar", 2, 1},
};

static const MolecularReactionSpecies D_2pAr2p_dissociation_products[] = {
    {"D", 1, 2},
    {"Ar", 1, 1},
};

static const MolecularReactionSpecies ArDp_dissociation_reactants[] = {
    {"e", -1, 1},
    {"ArD", 1, 1},
};

static const MolecularReactionSpecies ArDp_dissociation_products[] = {
    {"D", 0, 1},
    {"Ar", 0, 1},
};

static const MolecularReactionSpecies ArDp_dissociation_2_reactants[] = {
    {"e", -1, 1},
    {"ArD", 1, 1},
};

static const MolecularReactionSpecies ArDp_dissociation_2_products[] = {
    {"Ar", 1, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies ArDp_dissociation_3_reactants[] = {
    {"e", -1, 1},
    {"ArD", 1, 1},
};

static const MolecularReactionSpecies ArDp_dissociation_3_products[] = {
    {"Ar", 0, 1},
    {"D", 1, 1},
};

static const MolecularReactionSpecies Hep_He_charge_exchange_reactants[] = {
    {"He", 1, 1},
    {"He", 0, 1},
};

static const MolecularReactionSpecies Hep_He_charge_exchange_products[] = {
    {"He", 0, 1},
    {"He", 1, 1},
};

static const MolecularReactionSpecies D_Hep_charge_exchange_reactants[] = {
    {"D", 0, 1},
    {"He", 1, 1},
};

static const MolecularReactionSpecies D_Hep_charge_exchange_products[] = {
    {"D", 1, 1},
    {"He", 0, 1},
};

static const MolecularReactionSpecies Ar_Hep_charge_exchange_reactants[] = {
    {"Ar", 0, 1},
    {"He", 1, 1},
};

static const MolecularReactionSpecies Ar_Hep_charge_exchange_products[] = {
    {"Ar", 1, 1},
    {"He", 0, 1},
};

static const MolecularReactionSpecies HeDp_dissociation_reactants[] = {
    {"e", -1, 1},
    {"HeD", 1, 1},
};

static const MolecularReactionSpecies HeDp_dissociation_products[] = {
    {"He", 0, 1},
    {"D", 1, 1},
};

static const MolecularReactionSpecies He_D_2p_dissociation_reactants[] = {
    {"He", 0, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies He_D_2p_dissociation_products[] = {
    {"HeD", 1, 1},
    {"D", 0, 1},
};

static const MolecularReactionSpecies D_HeDp_dissociation_reactants[] = {
    {"D", 0, 1},
    {"HeD", 1, 1},
};

static const MolecularReactionSpecies D_HeDp_dissociation_products[] = {
    {"He", 0, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies D_2_HeDp_dissociation_reactants[] = {
    {"D2", 0, 1},
    {"HeD", 1, 1},
};

static const MolecularReactionSpecies D_2_HeDp_dissociation_products[] = {
    {"He", 0, 1},
    {"D3", 1, 1},
};

static const MolecularReactionSpecies Hep_D_2_charge_exchange_reactants[] = {
    {"He", 1, 1},
    {"D2", 1, 1},
};

static const MolecularReactionSpecies Hep_D_2_charge_exchange_products[] = {
    {"He", 0, 1},
    {"D2", 1, 1},
};

const MolecularReactionDefinition molecularReactionDefinitions[] = {
    {
        "Ar_Ar+_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Ar_Arp_charge_exchange_reactants,
        2, Ar_Arp_charge_exchange_products
    },
    {
        "Ar_Ar2+_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Ar_Ar2p_charge_exchange_reactants,
        2, Ar_Ar2p_charge_exchange_products
    },
    {
        "D_D+_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, D_Dp_charge_exchange_reactants,
        2, D_Dp_charge_exchange_products
    },
    {
        "D_2_D_2+_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, D_2_D_2p_charge_exchange_reactants,
        2, D_2_D_2p_charge_exchange_products
    },
    {
        "D_2+_D_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, D_2p_D_charge_exchange_reactants,
        2, D_2p_D_charge_exchange_products
    },
    {
        "Ar+_D_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Arp_D_charge_exchange_reactants,
        2, Arp_D_charge_exchange_products
    },
    {
        "Ar2+_D_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Ar2p_D_charge_exchange_reactants,
        2, Ar2p_D_charge_exchange_products
    },
    {
        "Ar3+_D_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Ar3p_D_charge_exchange_reactants,
        2, Ar3p_D_charge_exchange_products
    },
    {
        "Ar4+_D_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Ar4p_D_charge_exchange_reactants,
        2, Ar4p_D_charge_exchange_products
    },
    {
        "D_2_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, D_2_dissociation_reactants,
        2, D_2_dissociation_products
    },
    {
        "D_2_ionization",
        MolecularReactionProcess::IONIZATION,
        2, D_2_ionization_reactants,
        2, D_2_ionization_products
    },
    {
        "D_2_dissociation_2",
        MolecularReactionProcess::DISSOCIATIVE_IONIZATION,
        2, D_2_dissociation_2_reactants,
        2, D_2_dissociation_2_products
    },
    {
        "D+_D_2_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Dp_D_2_charge_exchange_reactants,
        2, Dp_D_2_charge_exchange_products
    },
    {
        "D_2+_D_2_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, D_2p_D_2_dissociation_reactants,
        2, D_2p_D_2_dissociation_products
    },
    {
        "Ar+_D_2_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Arp_D_2_charge_exchange_reactants,
        2, Arp_D_2_charge_exchange_products
    },
    {
        "Ar+_D_2_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, Arp_D_2_dissociation_reactants,
        2, Arp_D_2_dissociation_products
    },
    {
        "Ar2+_D_2_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Ar2p_D_2_charge_exchange_reactants,
        2, Ar2p_D_2_charge_exchange_products
    },
    {
        "Ar3+_D_2_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Ar3p_D_2_charge_exchange_reactants,
        2, Ar3p_D_2_charge_exchange_products
    },
    {
        "Ar4+_D_2_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Ar4p_D_2_charge_exchange_reactants,
        2, Ar4p_D_2_charge_exchange_products
    },
    {
        "ArD+_D_2_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, ArDp_D_2_dissociation_reactants,
        2, ArDp_D_2_dissociation_products
    },
    {
        "D_2+_Ar_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, D_2p_Ar_dissociation_reactants,
        2, D_2p_Ar_dissociation_products
    },
    {
        "D_3+_Ar_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, D_3p_Ar_dissociation_reactants,
        2, D_3p_Ar_dissociation_products
    },
    {
        "D_2+_recombination",
        MolecularReactionProcess::RECOMBINATION,
        2, D_2p_recombination_reactants,
        2, D_2p_recombination_products
    },
    {
        "D_2+_recombination_2",
        MolecularReactionProcess::RECOMBINATION,
        2, D_2p_recombination_2_reactants,
        2, D_2p_recombination_2_products
    },
    {
        "D_3+_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, D_3p_dissociation_reactants,
        2, D_3p_dissociation_products
    },
    {
        "D_3+_dissociation_2",
        MolecularReactionProcess::DISSOCIATION,
        2, D_3p_dissociation_2_reactants,
        1, D_3p_dissociation_2_products
    },
    {
        "D_2+Ar2+_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, D_2pAr2p_dissociation_reactants,
        2, D_2pAr2p_dissociation_products
    },
    {
        "ArD+_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, ArDp_dissociation_reactants,
        2, ArDp_dissociation_products
    },
    {
        "ArD+_dissociation_2",
        MolecularReactionProcess::DISSOCIATION,
        2, ArDp_dissociation_2_reactants,
        2, ArDp_dissociation_2_products
    },
    {
        "ArD+_dissociation_3",
        MolecularReactionProcess::DISSOCIATION,
        2, ArDp_dissociation_3_reactants,
        2, ArDp_dissociation_3_products
    },
    {
        "He+_He_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Hep_He_charge_exchange_reactants,
        2, Hep_He_charge_exchange_products
    },
    {
        "D_He+_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, D_Hep_charge_exchange_reactants,
        2, D_Hep_charge_exchange_products
    },
    {
        "Ar_He+_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Ar_Hep_charge_exchange_reactants,
        2, Ar_Hep_charge_exchange_products
    },
    {
        "HeD+_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, HeDp_dissociation_reactants,
        2, HeDp_dissociation_products
    },
    {
        "He_D_2+_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, He_D_2p_dissociation_reactants,
        2, He_D_2p_dissociation_products
    },
    {
        "D_HeD+_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, D_HeDp_dissociation_reactants,
        2, D_HeDp_dissociation_products
    },
    {
        "D_2_HeD+_dissociation",
        MolecularReactionProcess::DISSOCIATION,
        2, D_2_HeDp_dissociation_reactants,
        2, D_2_HeDp_dissociation_products
    },
    {
        "He+_D_2_charge_exchange",
        MolecularReactionProcess::CHARGE_EXCHANGE,
        2, Hep_D_2_charge_exchange_reactants,
        2, Hep_D_2_charge_exchange_products
    },
};

const len_t molecularReactionDefinitionCount =
    sizeof(molecularReactionDefinitions) / sizeof(molecularReactionDefinitions[0]);

} // namespace DREAM
