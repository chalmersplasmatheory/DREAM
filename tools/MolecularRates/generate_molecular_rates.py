
#!/usr/bin/env python3

import argparse
import datetime
import pathlib
import data
import numpy as np
from math import factorial
from scipy.special import gamma as gamma_function
from scipy.optimize import brentq



def cpp_array(realtype, name, values):
      s = f"const {realtype} {name}[{len(values)}] = {{"
      s += ",".join(f"{float(v):.12e}" for v in values)
      s += "};\n"
      return s


def integral_formula_eval(I_energy,m,l,v_ms,T,n, mass):
    ###this is still in cross section per velocity but we want to convert it to <sigma v> in m^3/s
     #I_energy is the ionization energy in eV
    Eh_ev = 27.211386245988  # Hartree energy in eV
    v0 = 2.1876912633e6  # Bohr velocity in m/s
    a0 = 5.29177210903e-11  # Bohr radius in m

    ##currently overriding v_ms with the thermal velocity, but we should use the actual relative 
    # velocity of the reactants
    E = 3/2 * T * n #remeber these are in log so change
    v_ms = (2*E/mass)**0.5 
    v_au = v_ms / v0

    I = I_energy/Eh_ev
    gamma = np.sqrt(2*I)
    nu = 1/gamma

    A = 1/np.sqrt(nu**2 *gamma_function(2*nu))

    def C_nulm(m, l):
        """ Calculate the C_{nu,l,m} coefficient. """
        return (
            (2*l + 1)
            * factorial(l + m)
            / (factorial(m) * factorial(l - m))
            * (-1)**(l + m + 1)
            / 2**(m + 1)
            * (4 / np.e)**nu
            * A**2
        )


    def Gamma_b(b,v_au,m, l):
         return (
            1 / v_au
            * np.sqrt(2*np.pi / (gamma*b))
            * C_nulm(m, l)
            * (gamma*b)**(2*nu - m)
            * np.exp(-gamma*b)
        )
    
    def find_b1(v_au, m, l):
        # Gamma(b) has a maximum at:
        #
        # b_peak = (2 nu - m - 1/2)/gamma
        #
        p = 2*nu - m - 0.5
        b_peak = p / gamma

        target = 1 / np.pi

        # Check that Gamma actually reaches 1/pi
        if abs(Gamma_b(b_peak, v_au, m, l)) <= target:
            raise ValueError(
                f"No b1 solution for m={m} at v={v_au:.3e} a.u."
            )

        def f(b):
            return abs(Gamma_b(b, v_au, m, l)) - target

        # Start at the peak and move outward until
        # Gamma has dropped below 1/pi.
        b_low = b_peak
        b_high = 2 * b_peak

        while f(b_high) > 0:
            b_high *= 2

        return brentq(f, b_low, b_high)
    def cross_section_m(v_au, m, l):

        b1 = find_b1(v_au, m, l)

        sigma_au = 0.5 * np.pi * b1**2

        return sigma_au

    #spefici for this
    sigma_m0 = cross_section_m(v_au, m, l)
    sigma_m1 = cross_section_m(v_au, l, l)

    sigma_au = (sigma_m0 + 2*sigma_m1) / 3

    # Convert a0^2 -> m^2
    sigma_tot = sigma_au * a0**2
    #Needs to be comvertet to <sigma v> in m^3/s correctly
    return sigma_tot * v_ms


    
         
         

def amjuel_eval(coeff, n, T, nT=9, nn=9):
      n_eval = min(max(n, 1e14), 1e22)
      logn = np.log(n_eval * 1e-14)   # n / 1e14 = n[m^-3] / (1e8 cm^-3)
      logT = np.log(T)

      s = 0.0
      pT = 1.0
      for iT in range(nT):
          pn = 1.0
          for in_ in range(nn):
              s += coeff[iT*nn + in_] * pT * pn
              pn *= logn
          pT *= logT

      return np.exp(s) * 1e-6         # cm^3/s -> m^3/s

def temperature_equation_eval(T):
      #turn T into Kelvin
      T = T * 11604.5250061657
     
      return (-1.3e-8 + 1.27e-6*T**(-0.48)) *1e-6  # cm^3/s -> m^3/s

def constant_eval(value, n, T):
      return value


def zero_eval(n, T):
      return 0.0

def build_rate_table(reaction, n_grid, T_grid):
      dtype = reaction["data_type"]
      raw = reaction["data"]

      values = []
      for T in T_grid: 
            for n in n_grid:
          
                if dtype in ("TODO", "FORMULA", "UNCLEAR"):
                    v = zero_eval(n, T)

                elif dtype == "CONSTANT":
                    v = constant_eval(raw, n, T)

                elif dtype == "INTEGRAL_FORMULA":
                    v = integral_formula_eval(
                        raw["I_energy"],
                        raw["m"],
                        raw["l"],
                        raw["v_ms"],
                        T,
                        n,
                        raw["mass"]
                    )
                elif dtype == "TEMPERATURE_EQUATION":
                    v = temperature_equation_eval(T)

                elif dtype in ("AMJUEL", "AMJUEL_POLYNOMIAL"):
                    coeff = getattr(data, raw["coefficients"])
                    v = amjuel_eval(
                        coeff,
                        n,
                        T,
                        raw.get("n_temperature", 9),
                        raw.get("n_density", 9)
                    )

                elif dtype == "ENERGY_INTERPOLATE":
                    # For now: placeholder zero until you decide the physics conversion.
                    # Later this becomes cross-section -> <sigma v>(n,T).
                    v = zero_eval(n, T)

                else:
                    raise ValueError(f"Unknown data_type '{dtype}' for {reaction['name']}")

                if v > 0:
                    values.append(np.log10(v))
                else:
                    values.append(-300)
                    #print(f"WARNING: Molecular rate '{reaction['name']}' evaluated to a negative value at n={n:.3e}, T={T:.3e}. Using -300.0 instead.")

      return values


def compile_molecular_rates(outputfile, inttype="len_t", realtype="real_t"):
      n_grid = np.logspace(14, 22, 50)   # m^-3
      T_grid = np.logspace(0, 4, 80)     # eV

      logn = np.log10(n_grid)
      logT = np.log10(T_grid)

      body = ""
      table = (f"const {inttype} molecular_rate_n = {len(data.MOLECULAR_REACTIONS)};\n"f"struct molecular_rate molecular_rate_table[{len(data.MOLECULAR_REACTIONS)}] ={{\n"
      )

      body += f"const {inttype} molecular_rate_nn = {len(n_grid)};\n"
      body += f"const {inttype} molecular_rate_nT = {len(T_grid)};\n"
      body += cpp_array(realtype, "molecular_rate_logn", logn)
      body += cpp_array(realtype, "molecular_rate_logT", logT)

      for reaction in data.MOLECULAR_REACTIONS:
          name = reaction["name"]
          cname = name.replace("+", "p").replace("-", "m").replace(".", "_")

          coeff_name = f"{cname}_coeff"
          coeff = build_rate_table(reaction, n_grid, T_grid)

          body += f"\n/* {name} */\n"
          body += cpp_array(realtype, coeff_name, coeff)

          table += (
              f'    {{"{name}", molecular_rate_nn, molecular_rate_nT, '
              f"molecular_rate_logn, molecular_rate_logT, {coeff_name}}},\n"
          )

      table += "};\n"

      filecontents = (
          "/* This file was auto-generated by 'generate_molecular_rates.py' "
          f"on {datetime.datetime.now().isoformat(sep=' ', timespec='seconds')} */\n\n"
          '#include "DREAM/MolecularRateData.hpp"\n\n'
          "namespace DREAM {\n\n"
          + body
          + "\n"
          + table
          + "\n}\n"
      )

      pathlib.Path(outputfile).parent.mkdir(parents=True, exist_ok=True)
      pathlib.Path(outputfile).write_text(filecontents)



def main():
      root = pathlib.Path(__file__).resolve().parents[2]
      outputfile = root / "src/MolecularRateData.cpp"
      compile_molecular_rates(outputfile)


if __name__ == "__main__":
      main()
