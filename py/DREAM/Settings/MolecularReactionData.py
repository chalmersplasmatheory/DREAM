from .MolecularReactionDataCoefficients import *


# Reactants/products: (name, charge, number of particles)
#
# Input selectors describe coefficient dependencies.
# kind="relative" with two species means their arithmetic-mean temperature
# Temperatures: eV; densities: m^-3.

# Charge Z of every species appearing in the reactions (summed nuclear
# charges for molecules).
SPECIES_Z = {
    "D":   1,
    "D2":  2,
    "D3":  3,
    "He":  2,
    "HeD": 3,
    "Ar":  18,
    "ArD": 19,
}

# Default mass [kg] of the molecules (Ions.addMolecule())
MOLECULE_MASS = {
    "D2":  2 * 3.3435837724e-27,
    "T2":  2 * 5.00826765e-27,
    "D3":  3 * 3.3435837724e-27,
    "ArD": 6.97e-26,
    "HeD": 9.99e-27,
    "DT":  8.35e-27,
}

# Delta E is computed from the potential energy E of each species, relative to
# neutral atoms:
#   D 0,  
#   D+ 13.60,  (NIST)
#   D2 -4.56,  
#   D2+ 10.91 (-4.56 + 15.47),
#   D3+ 4.64 (Hollmann 2019, eq. A.23d-f),
#   Ar 0,  
#   Ar+ 15.76,  
#   Ar2+ 43.39,  
#   Ar3+ 84.13,  
#   Ar4+ 143.71 (NIST),
#   ArD+ 9.70 (Hollmann 2019, eq. A.23i: Ar+ + D2 -> ArD+ + D releases 1.5),
#   He 0,  
#   He+ 24.59,  
#   HeD+ 11.76 (13.60 - 1.84, proton affinity of He).
#
# "franck_condon" [eV] (optional, default 0): kinetic energy that the fragments
# get when an electron breaks a molecule apart. The electron pays it on top of
# -Delta E (RadiatedPowerTerm), and the heavy products receive it
# (ReactionEnergyTransferTerm). Only given for the electron-impact breakups with
# a known value; e.g. D_2+_dissociation_2 and ArD+_dissociation_2/_3 use 0.
MOLECULAR_REACTIONS = [
    {
        "name": "Ar_Ar+_charge_exchange",
        "reactants": [("Ar", 0, 1), ("Ar", 1, 1)],
        "products": [("Ar", 1, 1), ("Ar", 0, 1)],
        "classification": "CHARGE_EXCHANGE_RESONANT",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar_Arp_cx_constant, "denominator": A_Ar},
        "temperature_input": {"kind": "species", "name": "Ar", "charge": 1},
        "density_input": {"kind": "none"},
        "energy_released": 0.0
    },

    {
        "name": "Ar_Ar2+_charge_exchange",
        "reactants": [("Ar", 0, 1), ("Ar", 2, 1)],
        "products": [("Ar", 2, 1), ("Ar", 0, 1)],
        "classification": "CHARGE_EXCHANGE_RESONANT",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar_Ar2p_cx_constant, "denominator": A_Ar},
        "temperature_input": {"kind": "species", "name": "Ar", "charge": 2},
        "density_input": {"kind": "none"},
        "energy_released": 0.0
    },

    {
        "name": "D_D+_charge_exchange",
        "reactants": [("D", 0, 1), ("D", 1, 1)],
        "products": [("D", 1, 1), ("D", 0, 1)],
        "classification": "CHARGE_EXCHANGE_RESONANT",
        "data_type": "APPROXIMATION",
        "data": {"constant": D_Dp_cx_constant, "denominator": A_D},
        "temperature_input": {"kind": "species", "name": "D", "charge": 1},
        "density_input": {"kind": "none"},
        "energy_released": 0.0
    },

    {
        "name": "D_2_D_2+_charge_exchange",
        "reactants": [("D2", 0, 1), ("D2", 1, 1)],
        "products": [("D2", 1, 1), ("D2", 0, 1)],
        "classification": "CHARGE_EXCHANGE_RESONANT",
        "data_type": "APPROXIMATION",
        "data": {"constant": D2_D2p_cx_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 1},
        "density_input": {"kind": "none"},
        "energy_released": 0.0
    },

    {
        "name": "D_2+_D_charge_exchange",
        "reactants": [("D2", 1, 1), ("D", 0, 1)],
        "products": [("D2", 0, 1), ("D", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "CONSTANT",
        "data": D2p_D_cx_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": 1.87
    },

    {
        "name": "Ar+_D_charge_exchange",
        "reactants": [("Ar", 1, 1), ("D", 0, 1)],
        "products": [("Ar", 0, 1), ("D", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar_D_cx_constant_1, "denominator": A_D},
        "temperature_input": {"kind": "species", "name": "D", "charge": 0},
        "density_input": {"kind": "none"},
        "energy_released": 2.16
    },

    {
        "name": "Ar2+_D_charge_exchange",
        "reactants": [("Ar", 2, 1), ("D", 0, 1)],
        "products": [("Ar", 1, 1), ("D", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar_D_cx_constant_2, "denominator": A_D},
        "temperature_input": {"kind": "species", "name": "D", "charge": 0},
        "density_input": {"kind": "none"},
        "energy_released": 14.03
    },

    {
        "name": "Ar3+_D_charge_exchange",
        "reactants": [("Ar", 3, 1), ("D", 0, 1)],
        "products": [("Ar", 2, 1), ("D", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar_D_cx_constant_3, "denominator": A_D},
        "temperature_input": {"kind": "species", "name": "D", "charge": 0},
        "density_input": {"kind": "none"},
        "energy_released": 27.14
    },

    {
        "name": "Ar4+_D_charge_exchange",
        "reactants": [("Ar", 4, 1), ("D", 0, 1)],
        "products": [("Ar", 3, 1), ("D", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar_D_cx_constant_4, "denominator": A_D},
        "temperature_input": {"kind": "species", "name": "D", "charge": 0},
        "density_input": {"kind": "none"},
        "energy_released": 45.98
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
        "density_input": {"kind": "electron"},
        "energy_released": -4.56,
        "franck_condon": 5.38   # AMJUEL H.10 2.2.h2c: 9.94 - 4.56 (Hollmann 2019, eq. A.23a: 5)
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
        "density_input": {"kind": "electron"},
        "energy_released": -15.47
    },

    {
        "name": "D_2_dissociation_2",
        "reactants": [("e", -1, 1), ("D2", 0, 1)],
        "products": [("e", -1, 2), ("D", 0, 1), ("D", 1, 1)],
        "classification": "DISSOCIATION",
        "data_type": "AMJUEL_POLYNOMIAL",
        "data": {
            "n_temperature": 9,
            "n_density": 9,
            "coefficients": "D2_dissioniz_2210_coeff"
        },
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"},
        "energy_released": -18.16,
        "franck_condon": 6.0   # Hollmann 2019, eq. A.23b
    },

    {
        "name": "D_2_runaway",
        "reactants": [("RE", -1, 1), ("D2", 0, 1)],
        "products": [("RE", -1, 1), ("D2", 1, 1), ("e", -1, 1)],
        "classification": "RUNAWAY",
        "data_type": "RIEKE",
        "data": {"M2": D2_rieke_M2, "C": D2_rieke_C, "factor": 1},
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": -15.47
    },

    {
        "name": "D+_D_2_charge_exchange",
        "reactants": [("D", 1, 1), ("D2", 0, 1)],
        "products": [("D", 0, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "UNIQUE_FORMULA_1",
        "data": {},
        "temperature_input": {"kind": "species", "name": "D", "charge": 1},
        "density_input": {"kind": "none"},
        "energy_released": -1.87
    },

    {
        "name": "D_2+_D_2_dissociation",
        "reactants": [("D2", 1, 1), ("D2", 0, 1)],
        "products": [("D3", 1, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": D2p_D2_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": 1.71
    },

    {
        "name": "Ar+_D_2_charge_exchange",
        "reactants": [("Ar", 1, 1), ("D2", 0, 1)],
        "products": [("Ar", 0, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Arp_D2_cx_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 0},
        "density_input": {"kind": "none"},
        "energy_released": 0.29
    },

    {
        "name": "Ar+_D_2_dissociation",
        "reactants": [("Ar", 1, 1), ("D2", 0, 1)],
        "products": [("ArD", 1, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": Arp_D2_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": 1.5
    },

    {
        "name": "Ar2+_D_2_charge_exchange",
        "reactants": [("Ar", 2, 1), ("D2", 0, 1)],
        "products": [("Ar", 1, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar2p_D2_cx_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 0},
        "density_input": {"kind": "none"},
        "energy_released": 12.16
    },

    {
        "name": "Ar3+_D_2_charge_exchange",
        "reactants": [("Ar", 3, 1), ("D2", 0, 1)],
        "products": [("Ar", 2, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar3p_D2_cx_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 0},
        "density_input": {"kind": "none"},
        "energy_released": 25.27
    },

    {
        "name": "Ar4+_D_2_charge_exchange",
        "reactants": [("Ar", 4, 1), ("D2", 0, 1)],
        "products": [("Ar", 3, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "APPROXIMATION",
        "data": {"constant": Ar4p_D2_cx_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 0},
        "density_input": {"kind": "none"},
        "energy_released": 44.11
    },

    {
        "name": "ArD+_D_2_dissociation",
        "reactants": [("ArD", 1, 1), ("D2", 0, 1)],
        "products": [("Ar", 0, 1), ("D3", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": ArDp_D2_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": 0.5
    },

    {
        "name": "D_2+_Ar_dissociation",
        "reactants": [("D2", 1, 1), ("Ar", 0, 1)],
        "products": [("D", 0, 1), ("ArD", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": D2p_Ar_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": 1.21
    },

    {
        "name": "D_3+_Ar_dissociation",
        "reactants": [("D3", 1, 1), ("Ar", 0, 1)],
        "products": [("D2", 0, 1), ("ArD", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": D3p_Ar_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": -0.5
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
        "density_input": {"kind": "electron"},
        "energy_released": 10.91
    },

    {
        "name": "D_2+_dissociation",
        "reactants": [("e", -1, 1), ("D2", 1, 1)],
        "products": [("e", -1, 1), ("D", 0, 1), ("D", 1, 1)],
        "classification": "DISSOCIATION",
        "data_type": "AMJUEL_POLYNOMIAL",
        "data": {
            "n_temperature": 9,
            "n_density": 9,
            "coefficients": "D2_diss_2212_coeff"
        },
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"},
        "energy_released": -2.69,
        "franck_condon": 7.0   # Hollmann 2019, eq. A.23c
    },

    {
        "name": "D_2+_dissociation_2",
        "reactants": [("e", -1, 1), ("D2", 1, 1)],
        "products": [("e", -1, 2), ("D", 1, 2)],
        "classification": "DISSOCIATION",
        "data_type": "AMJUEL_POLYNOMIAL",
        "data": {
            "n_temperature": 9,
            "n_density": 9,
            "coefficients": "D2_dissioniz_2211_coeff"
        },
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "electron"},
        "energy_released": -16.29
    },

    {
        "name": "D_2+_runaway",
        "reactants": [("RE", -1, 1), ("D2", 1, 1)],
        "products": [("RE", -1, 1), ("D", 1, 2), ("e", -1, 1)],
        "classification": "RUNAWAY",
        "data_type": "RIEKE",
        "data": {"M2": H_rieke_M2, "C": H_rieke_C, "factor": 1},
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": -16.29
    },

    {
        "name": "D_3+_dissociation",
        "reactants": [("e", -1, 1), ("D3", 1, 1)],
        "products": [("D2", 0, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION",
        "data_type": "TEMPERATURE_EQUATION",
        "data": D3p_diss_D2_D_fraction,
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"},
        "energy_released": 9.2
    },

    {
        "name": "D_3+_dissociation_2",
        "reactants": [("e", -1, 1), ("D3", 1, 1)],
        "products": [("D", 0, 3)],
        "classification": "DISSOCIATION",
        "data_type": "TEMPERATURE_EQUATION",
        "data": D3p_diss_3D_fraction,
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"},
        "energy_released": 4.64
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
        "density_input": {"kind": "none"},
        "energy_released": -19.87
    },

    {
        "name": "D_2+Ar2+_dissociation",
        "reactants": [("D2", 1, 1), ("Ar", 2, 1)],
        "products": [("D", 1, 2), ("Ar", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "APPROXIMATION",
        "data": {"constant": D2p_Ar2p_diss_constant, "denominator": A_D2},
        "temperature_input": {"kind": "species", "name": "D2", "charge": 1},
        "density_input": {"kind": "none"},
        "energy_released": 11.34
    },

    {
        "name": "ArD+_dissociation",
        "reactants": [("e", -1, 1), ("ArD", 1, 1)],
        "products": [("D", 0, 1), ("Ar", 0, 1)],
        "classification": "DISSOCIATION",
        "data_type": "UNIQUE_FORMULA_2",
        "data": {},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"},
        "energy_released": 9.7
    },

    {
        "name": "ArD+_dissociation_2",
        "reactants": [("e", -1, 1), ("ArD", 1, 1)],
        "products": [("e", -1, 1), ("Ar", 1, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION",
        "data_type": "UNIQUE_FORMULA_3",
        "data": {},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"},
        "energy_released": -6.06
    },

    {
        "name": "ArD+_dissociation_3",
        "reactants": [("e", -1, 1), ("ArD", 1, 1)],
        "products": [("e", -1, 1), ("Ar", 0, 1), ("D", 1, 1)],
        "classification": "DISSOCIATION",
        "data_type": "UNIQUE_FORMULA_3",
        "data": {},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"},
        "energy_released": -3.9
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
        "density_input": {"kind": "none"},
        "energy_released": -19.66
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
        "density_input": {"kind": "none"},
        "energy_released": 0.0
    },

    {
        "name": "D_He+_charge_exchange",
        "reactants": [("D", 0, 1), ("He", 1, 1)],
        "products": [("D", 1, 1), ("He", 0, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "POLYFIT_2",
        "data": D_He_cx,
        "temperature_input": {"kind": "species", "name": "He", "charge": 1},
        "density_input": {"kind": "none"},
        "energy_released": 10.99
    },

    {
        "name": "Ar_He+_charge_exchange",
        "reactants": [("Ar", 0, 1), ("He", 1, 1)],
        "products": [("Ar", 1, 1), ("He", 0, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "UNIQUE_FORMULA_4",
        "data": {},
        "temperature_input": {"kind": "species", "name": "He", "charge": 1},
        "density_input": {"kind": "none"},
        "energy_released": 8.83
    },

    {
        "name": "HeD+_dissociation",
        "reactants": [("e", -1, 1), ("HeD", 1, 1)],
        "products": [("He", 0, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "UNIQUE_FORMULA_5",
        "data": {},
        "temperature_input": {"kind": "electron"},
        "density_input": {"kind": "none"},
        "energy_released": 11.76
    },

    {
        "name": "He_D_2+_dissociation",
        "reactants": [("He", 0, 1), ("D2", 1, 1)],
        "products": [("HeD", 1, 1), ("D", 0, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": He_D2p_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": -0.85
    },

    {
        "name": "D_HeD+_dissociation",
        "reactants": [("D", 0, 1), ("HeD", 1, 1)],
        "products": [("He", 0, 1), ("D2", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": D_HeDp_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": 0.85
    },

    {
        "name": "D_2_HeD+_dissociation",
        "reactants": [("D2", 0, 1), ("HeD", 1, 1)],
        "products": [("He", 0, 1), ("D3", 1, 1)],
        "classification": "DISSOCIATION_MULTIPLE",
        "data_type": "CONSTANT",
        "data": D2_HeDp_diss_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": 2.56
    },

    {
        "name": "He+_D_2_charge_exchange",
        "reactants": [("He", 1, 1), ("D2", 0, 1)],
        "products": [("He", 0, 1), ("D2", 1, 1)],
        "classification": "CHARGE_EXCHANGE",
        "data_type": "CONSTANT",
        "data": Hep_D2_cx_rate,
        "temperature_input": {"kind": "none"},
        "density_input": {"kind": "none"},
        "energy_released": 9.12
    },
]



# Line radiation of molecular species: the molecular counterpart of ADAS PLT.
# Not reactions (no density changes); used in RadiatedPowerTerm as
#   P_rad = n_e * n(species, Z0) * L(n_e, T_e),   L in eV m^3/s.
# Used whenever the species is in the simulation and molecular reactions are enabled.
MOLECULAR_LINE_RADIATION = [
    {
        # AMJUEL H.10 2.2.h2r: radiation from e + D2 (dissociation + ionization)
        "name": "D_2_radiation",
        "species": ("D2", 0),
        "data_type": "AMJUEL_POLYNOMIAL",
        "data": {
            "n_temperature": 9,
            "n_density": 9,
            "coefficients": "D2_radiation_h2r_coeff"
        }
    },
]
