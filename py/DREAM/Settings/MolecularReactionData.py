from .data_coeff import *


# Reactants/products: (name, charge, number of particles)
#
# Input selectors describe coefficient depandancies.
# kind="relative" with two species means their arithmetic-mean temperature
# Temperatures: eV; densities: m^-3.
MOLECULAR_REACTIONS = [
    {
        "name": "Ar_Ar+_charge_exchange",
        "reactants": [("Ar", 0, 1), ("Ar", 1, 1)],
        "products": [("Ar", 1, 1), ("Ar", 0, 1)],
        "classification": "CHARGE_EXCHANGE_RESONANT",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar_Arp_cx_constant, "denominator": A_Ar},
        "temperature_input": {"kind": "species", "name": "Ar", "charge": 1},
        "density_input": {"kind": "none"}
    },

    {
        "name": "Ar_Ar2+_charge_exchange",
        "reactants": [("Ar", 0, 1), ("Ar", 2, 1)],
        "products": [("Ar", 2, 1), ("Ar", 0, 1)],
        "classification": "CHARGE_EXCHANGE_RESONANT",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar_Ar2p_cx_constant, "denominator": A_Ar},
        "temperature_input": {"kind": "species", "name": "Ar", "charge": 2},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_D+_charge_exchange",
        "reactants": [("D", 0, 1), ("D", 1, 1)],
        "products": [("D", 1, 1), ("D", 0, 1)],
        "classification": "CHARGE_EXCHANGE_RESONANT",
        "data_type": "APPROXIMATION",
        "data": {"constant": D_Dp_cx_constant, "denominator": A_D},
        "temperature_input": {"kind": "species", "name": "D", "charge": 1},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_2_D_2+_charge_exchange",
        "reactants": [("D2", 0, 1), ("D2", 1, 1)],
        "products": [("D2", 1, 1), ("D2", 0, 1)],
        "classification": "CHARGE_EXCHANGE_RESONANT",
        "data_type": "APPROXIMATION",
        "data": {"constant": D2_D2p_cx_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 1},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_2+_D_charge_exchange",
        "reactants": [("D2", 1, 1), ("D", 0, 1)],
        "products": [("D2", 0, 1), ("D", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "CONSTANT",
        "data": D2p_D_cx_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "Ar+_D_charge_exchange",
        "reactants": [("Ar", 1, 1), ("D", 0, 1)],
        "products": [("Ar", 0, 1), ("D", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "ADAS",
        "data": {"filename": Ar_D_cx_adas_file, "Z0": 1},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"}
    },

    {
        "name": "Ar2+_D_charge_exchange",
        "reactants": [("Ar", 2, 1), ("D", 0, 1)],
        "products": [("Ar", 1, 1), ("D", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "ADAS",
        "data": {"filename": Ar_D_cx_adas_file, "Z0": 2},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"}
    },

    {
        "name": "Ar3+_D_charge_exchange",
        "reactants": [("Ar", 3, 1), ("D", 0, 1)],
        "products": [("Ar", 2, 1), ("D", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "ADAS",
        "data": {"filename": Ar_D_cx_adas_file, "Z0": 3},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"}
    },

    {
        "name": "Ar4+_D_charge_exchange",
        "reactants": [("Ar", 4, 1), ("D", 0, 1)],
        "products": [("Ar", 3, 1), ("D", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "ADAS",
        "data": {"filename": Ar_D_cx_adas_file, "Z0": 4},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"}
    },

    {
        "name": "D_2_dissociation",
        "reactants": [("e", -1, 1), ("D2", 0, 1)],
        "products": [("e", -1, 1), ("D", 0, 2)],
        "classification": "DISSOCIATION",
        "data_type": "AMJUEL_POLYNOMIAL",
        "data": {
            "n_temperature": 9,
            "n_density": 9,
            "coefficients": "D2_diss_225_coeff"
        },
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"}
    },

    {
        "name": "D_2_ionization",
        "reactants": [("e", -1, 1), ("D2", 0, 1)],
        "products": [("e", -1, 2), ("D2", 1, 1)],
        "classification": "IONIZATION",
        "data_type": "AMJUEL_POLYNOMIAL",
        "data": {
            "n_temperature": 9,
            "n_density": 9,
            "coefficients": "D2_ioniz_229_coeff"
        },
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"}
    },

    {
        "name": "D_2_dissociation_2",
        "reactants": [("e", -1, 1), ("D2", 0, 1)],
        "products": [("D", 0, 1), ("D", 1, 1)],
        "classification": "DISSOCIATION",
        "data_type": "AMJUEL_POLYNOMIAL",
        "data": {
            "n_temperature": 9,
            "n_density": 9,
            "coefficients": "D2_dissioniz_2210_coeff"
        },
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"}
    },

    {
        "name": "D_2_runaway",
        "reactants": [("RE", -1, 1), ("D2", 0, 1)],
        "products": [("RE", -1, 1), ("D2", 1, 1), ("e", -1, 1)],
        "classification": "RUNAWAY",
        "data_type": "RIEKE",
        "data": {"M2": D2_rieke_M2, "C": D2_rieke_C, "factor": 1},
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D+_D_2_charge_exchange",
        "reactants": [("D", 1, 1), ("D2", 0, 1)],
        "products": [("D", 0, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "UNIQUE_FORMULA_1",
        "data": {},
        "temperature_input": {"kind": "species", "name": "D", "charge": 1},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_2+_D_2_dissociation",
        "reactants": [("D2", 1, 1), ("D2", 0, 1)],
        "products": [("D3", 1, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": D2p_D2_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "Ar+_D_2_charge_exchange",
        "reactants": [("Ar", 1, 1), ("D2", 0, 1)],
        "products": [("Ar", 0, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Arp_D2_cx_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 0},
        "density_input": {"kind": "none"}
    },

    {
        "name": "Ar+_D_2_dissociation",
        "reactants": [("Ar", 1, 1), ("D2", 0, 1)],
        "products": [("ArD", 1, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": Arp_D2_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "Ar2+_D_2_charge_exchange",
        "reactants": [("Ar", 2, 1), ("D2", 0, 1)],
        "products": [("Ar", 1, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar2p_D2_cx_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 0},
        "density_input": {"kind": "none"}
    },

    {
        "name": "Ar3+_D_2_charge_exchange",
        "reactants": [("Ar", 3, 1), ("D2", 0, 1)],
        "products": [("Ar", 2, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar3p_D2_cx_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 0},
        "density_input": {"kind": "none"}
    },

    {
        "name": "Ar4+_D_2_charge_exchange",
        "reactants": [("Ar", 4, 1), ("D2", 0, 1)],
        "products": [("Ar", 3, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar4p_D2_cx_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 0},
        "density_input": {"kind": "none"}
    },

    {
        "name": "ArD+_D_2_dissociation",
        "reactants": [("ArD", 1, 1), ("D2", 0, 1)],
        "products": [("Ar", 0, 1), ("D3", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": ArDp_D2_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_2+_Ar_dissociation",
        "reactants": [("D2", 1, 1), ("Ar", 0, 1)],
        "products": [("D", 0, 1), ("ArD", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": D2p_Ar_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_3+_Ar_dissociation",
        "reactants": [("D3", 1, 1), ("Ar", 0, 1)],
        "products": [("D2", 0, 1), ("ArD", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": D3p_Ar_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_2+_recombination",
        "reactants": [("e", -1, 1), ("D2", 1, 1)],
        "products": [("D", 0, 1), ("D", 0, 1)],
        "classification": "RECOMBINATION",
        "data_type": "AMJUEL_POLYNOMIAL",
        "data": {
            "n_temperature": 9,
            "n_density": 9,
            "coefficients": "D2_diss_2214_coeff"
        },
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"}
    },

    {
        "name": "D_2+_dissociation",
        "reactants": [("e", -1, 1), ("D2", 1, 1)],
        "products": [("D", 0, 1), ("D", 1, 1)],
        "classification": "DISSOCIATION",
        "data_type": "AMJUEL_POLYNOMIAL",
        "data": {
            "n_temperature": 9,
            "n_density": 9,
            "coefficients": "D2_diss_2212_coeff"
        },
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"}
    },

    {
        "name": "D_2+_runaway",
        "reactants": [("RE", -1, 1), ("D2", 1, 1)],
        "products": [("RE", -1, 1), ("D", 1, 2), ("e", -1, 1)],
        "classification": "RUNAWAY",
        "data_type": "RIEKE",
        "data": {"M2": H_rieke_M2, "C": H_rieke_C, "factor": 1},
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_3+_dissociation",
        "reactants": [("e", -1, 1), ("D3", 1, 1)],
        "products": [("D2", 0, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION",
        "data_type": "TEMPERATURE_EQUATION",
        "data": D3p_diss_D2_D_fraction,
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_3+_dissociation_2",
        "reactants": [("e", -1, 1), ("D3", 1, 1)],
        "products": [("D", 0, 3)],
        "classification": "DISSOCIATION",
        "data_type": "TEMPERATURE_EQUATION",
        "data": D3p_diss_3D_fraction,
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_3+_runaway",
        "reactants": [("RE", -1, 1), ("D3", 1, 1)],
        "products": [("RE", -1, 1), ("D", 1, 1), ("D2", 1, 1), ("e", -1, 1)],
        "classification": "RUNAWAY",
        # 2 x H atom used as proxy (as in IonBalance)
        "data_type": "RIEKE",
        "data": {"M2": H_rieke_M2, "C": H_rieke_C, "factor": 2},
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_2+Ar2+_dissociation",
        "reactants": [("D2", 1, 1), ("Ar", 2, 1)],
        "products": [("D", 1, 2), ("Ar", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "APPROXIMATION",
        "data": {"constant": D2p_Ar2p_diss_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 1},
        "density_input": {"kind": "none"}
    },

    {
        "name": "ArD+_dissociation",
        "reactants": [("e", -1, 1), ("ArD", 1, 1)],
        "products": [("D", 0, 1), ("Ar", 0, 1)],
        "classification": "DISSOCIATION",
        "data_type": "UNIQUE_FORMULA_2",
        "data": {},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "ArD+_dissociation_2",
        "reactants": [("e", -1, 1), ("ArD", 1, 1)],
        "products": [("Ar", 1, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION",
        "data_type": "UNIQUE_FORMULA_3",
        "data": {},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "ArD+_dissociation_3",
        "reactants": [("e", -1, 1), ("ArD", 1, 1)],
        "products": [("Ar", 0, 1), ("D", 1, 1)],
        "classification": "DISSOCIATION",
        "data_type": "UNIQUE_FORMULA_3",
        "data": {},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "ArD+_runaway",
        "reactants": [("RE", -1, 1), ("ArD", 1, 1)],
        "products": [("RE", -1, 1), ("Ar", 1, 1), ("D", 1, 1), ("e", -1, 1)],
        "classification": "RUNAWAY",
        # Neutral Ar used as proxy (as in IonBalance)
        "data_type": "RIEKE",
        "data": {"M2": Ar_rieke_M2, "C": Ar_rieke_C, "factor": 1},
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "He+_He_charge_exchange",
        "reactants": [("He", 1, 1), ("He", 0, 1)],
        "products": [("He", 0, 1), ("He", 1, 1)],
        "classification": "CHARGE_EXCHANGE_RESONANT",
        "data_type": "POLYFIT",
        "data": He_cx,
        "temperature_input": {
            "kind": "relative",
            "species": [{"name": "He", "charge": 0}, {"name": "He", "charge": 1}]
        },
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_He+_charge_exchange",
        "reactants": [("D", 0, 1), ("He", 1, 1)],
        "products": [("D", 1, 1), ("He", 0, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "POLYFIT_2",
        "data": D_He_cx,
        "temperature_input": {"kind": "species", "name": "He", "charge": 1},
        "density_input": {"kind": "none"}
    },

    {
        "name": "Ar_He+_charge_exchange",
        "reactants": [("Ar", 0, 1), ("He", 1, 1)],
        "products": [("Ar", 1, 1), ("He", 0, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "UNIQUE_FORMULA_4",
        "data": {},
        "temperature_input": {"kind": "species", "name": "He", "charge": 1},
        "density_input": {"kind": "none"}
    },

    {
        "name": "HeD+_dissociation",
        "reactants": [("e", -1, 1), ("HeD", 1, 1)],
        "products": [("He", 0, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "UNIQUE_FORMULA_5",
        "data": {},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "He_D_2+_dissociation",
        "reactants": [("He", 0, 1), ("D2", 1, 1)],
        "products": [("HeD", 1, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": He_D2p_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_HeD+_dissociation",
        "reactants": [("D", 0, 1), ("HeD", 1, 1)],
        "products": [("He", 0, 1), ("D2", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": D_HeDp_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "D_2_HeD+_dissociation",
        "reactants": [("D2", 0, 1), ("HeD", 1, 1)],
        "products": [("He", 0, 1), ("D3", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": D2_HeDp_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },

    {
        "name": "He+_D_2_charge_exchange",
        "reactants": [("He", 1, 1), ("D2", 0, 1)],
        "products": [("He", 0, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "CONSTANT",
        "data": Hep_D2_cx_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"}
    },
]


