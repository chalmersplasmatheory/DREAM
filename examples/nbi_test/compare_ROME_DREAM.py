##Implemnteded from ROME paper: 
# "NEUTRAL-BEAM INJECTION INTO A TOKAMAK: Part I:
# FAST-ION SPATIAL DISTRIBUTION FOR TANGENTIAL INJECTION5"
# and then compared with DREAM results

#!/usr/bin/env python3
import h5py
import matplotlib.pyplot as plt
import numpy as np
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import trapezoid
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle


from DREAM.DREAMOutput import DREAMOutput

# Load simulation results
#do = DREAMOutput("output_nbi_50.h5")
do = DREAMOutput("output_nbi.h5")

# --- Plot T_cold as function of radius at final time ---
r = do.grid.r
# Constants and geometry for ROME
a = 0.23             # Minor radius [m]
R0 = 0.798            # Major radius [m]
xs = 0.0       # Stagnation point offset [m]
rB = 0.0775         # Beam radius [m]
Rc = 0.685  
lambda0= 0.434 #Taken from dream to match
jB = 250e3
from scipy.integrate import quad


def currentdensity(R_B):
    # Constant current density of 250 kA/m² within beam radius
    return jB

def calculate_IB():
    R_min = Rc - rB
    R_max = Rc + rB
    integrand = lambda R: currentdensity(R) * 2 * np.pi * R
    IB, _ = quad(integrand, R_min, R_max)

    return jB * np.pi * rB**2  

def U_func(x):
    return 0 if x < 0 else 1 if x > 0 else 0.5

def F_func(R, R_B, zB):
    value = np.sqrt(R**2-R_B**2)/a**2 * (a**2 -zB**2-R0**2 + R0*R -2/3*R_B**2-1/3*R**2) + R0*R_B**2/a**2 * np.log(R+np.sqrt(R**2-R_B**2))
    return value

def dNfdR_func_correct(R, R_B, zB):
    term3 = U_func(R - R_B)
    if term3 == 0 or R - R_B <= 1e-6:
        return 0
    f1=F_func(R, R_B, zB)
    f2 = F_func(R0 + np.sqrt(a**2-zB**2), R_B, zB)
    f3=F_func(R_B, R_B, zB)
    term2 = np.exp(1/lambda0)*(f2-f1) + np.exp(1/lambda0)*(f2-2*f3 + f1)
    term1 = ((a**2 - (R - R0)**2 - zB**2) / a**2) * (R / np.sqrt(R**2 - R_B**2)) * 1/lambda0

    return term1 * term2 * term3

def H_r(rho):
    if rho <= 0:
        return 0
    integral_1 = 0
    dz = min(rho, rB)/50
    integration_steps_1 = np.linspace(0, min(rho, rB), 50)
    for zB in integration_steps_1:
        if rho**2 - zB**2 <= 0:
            continue
        sqrt_rho = np.sqrt(rho**2 - zB**2)
        sqrt_rB = np.sqrt(rB**2 - zB**2)

        # Incoming beam half
        integral2_1 = 0
        R = R0 + xs + sqrt_rho
        R_B_min = Rc - sqrt_rB
        R_B_max = min(R, Rc + sqrt_rB)
        dR = (R_B_max - R_B_min)/50
        integration_steps_2_1 = np.linspace(R_B_min, R_B_max, 50)
        for R_B in integration_steps_2_1:
            integral2_1+=currentdensity(R_B) * dNfdR_func_correct(R, R_B, zB)*dR
      

        # Outgoing beam half
        integral2_2 = 0
        R = R0 + xs - sqrt_rho
        R_B_max = min(R, Rc + sqrt_rB)
        dR = (R_B_max - R_B_min)/50
        
        integration_steps_2_2 = np.linspace(R_B_min, R_B_max, 50)
        for R_B in integration_steps_2_2:
            integral2_2+=currentdensity(R_B) * dNfdR_func_correct(R, R_B, zB) *dR

        integral_1 += (integral2_1 + integral2_2) / sqrt_rho *dz
    return (a**2 / calculate_IB()) * integral_1/2

def H_r_theta(r, theta):
    rho = np.sqrt(r**2 + xs**2 - 2*r*xs*np.cos(theta)) 
    return H_r(rho) 

def real_H_r(r):
    theta_vals = np.linspace(0, 2*np.pi, 1)  # Does not matter when xs=0
    temp=[]
    for theta in theta_vals:
        temp.append(H_r_theta(r, theta))
    return np.mean(temp) /(2*np.pi)

def total_deposition_flux_surface(r):
    h_r = real_H_r(r)
    return h_r * 4 * np.pi**2 * R0 * r * a/len(dream_r)   


####################PLOTTING######################
dream_r = do.grid.r

rome_total = np.array([real_H_r(r) for r in dream_r])
rome_fraction_sum = np.array([total_deposition_flux_surface(r) for r in dream_r])
print("Fraction of deposition ROME:", np.sum([total_deposition_flux_surface(r) for r in dream_r])) 

H_r = do.other.fluid.Tcold_NBI[-1, :]


#TOTAL DEPOSITION PER FLUX SURFACE
plt.plot(dream_r[1:len(dream_r)], rome_fraction_sum[1:len(dream_r)]/(a / len(dream_r)), label='ROME', linewidth=2)
H_r_deposition=[]
for i in range(len(dream_r)):
    H_r_deposition.append(-H_r[i] * 4 * np.pi**2 * R0 * dream_r[i])
plt.scatter(do.grid.r[1:len(dream_r)], H_r_deposition[1:len(dream_r)], label = "DREAM", color='#00cc44') 
plt.xlabel(r'$r$ [m]') 
plt.ylabel(r'$H(r) V^\prime [1/$m$]$')
plt.legend()
plt.savefig("ROMECOMPARISON_H_depostionperflux.pdf", bbox_inches="tight")
plt.show()


#TOTAL FRACTION PER FLUX SURFACE
H_r_fraction = []
plt.plot(dream_r[1:len(dream_r)], rome_fraction_sum[1:len(dream_r)], label='ROME',linewidth=2)
for i in range(len(dream_r)):
    H_r_fraction.append(-H_r[i] * 4 * np.pi**2 * R0 * a / len(dream_r) * dream_r[i])
plt.scatter(do.grid.r[1:len(dream_r)], H_r_fraction[1:len(dream_r)], label = "DREAM", color='#00cc44') 

plt.xlabel(r'$r$ [m]') 
plt.ylabel(r'$H(r) V^\prime$ [1/$m$]')
plt.legend()
plt.savefig("ROMECOMPARISON_H_fraction.pdf", bbox_inches="tight")
plt.show()

#H(r) [m-3]
rome_H_r = []
dream_H_r = []
for i in range(len(dream_r)):
    dream_H_r.append(-H_r[i])

for i in range(len(dream_r)):
    rome_H_r.append(rome_total[i])

plt.plot(dream_r[1:len(dream_r)], rome_H_r[1:len(dream_r)], label='ROME', linewidth=2)
plt.scatter(do.grid.r[1:len(dream_r)], dream_H_r[1:len(dream_r)], label='DREAM', color='#00cc44')
plt.xlabel(r'$r$ [m]')
plt.ylabel(r'$H(r)$[1/m$^3$]')
plt.legend()
plt.savefig("ROMECOMPARISON_H.pdf", bbox_inches="tight")
plt.show()


##CUMULATED DEPOSITION
plt.plot(dream_r[1:len(dream_r)], np.cumsum(rome_fraction_sum[1:len(dream_r)]), label='ROME', linewidth=2)
plt.scatter(do.grid.r[1:len(dream_r)], np.cumsum(H_r_fraction[1:len(dream_r)]), label='DREAM', color='#00cc44')
plt.xlabel(r'$r$ [m]')
plt.ylabel(r'$\int_0^r H(r) V^\prime dr$')
plt.legend()
plt.savefig("ROMECOMPARISON_CUMULATED.pdf", bbox_inches="tight")
plt.show()
