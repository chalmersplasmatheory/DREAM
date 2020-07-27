#!/usr/bin/env python3
#
# This example first sets up a simple runaway scenario, which is
# then passed to a DREAM.ConvergenceScan object. The convergence
# scan is configured to apply to the most relevant resolution
# parameters for this scenario. We will use the runaway rate as
# a measure of convergence, and we will consider it to be converged
# when it no longer varies significantly. Once set up, the convergence
# scan is also executed and its result are presented.
#
# Run as
#
#   $ ./run.py
#
# ###################################################################

import numpy as np
import sys

sys.path.append('../../py/')

from DREAM.ConvergenceScanPlot import ConvergenceScanPlot
from DREAM import DREAMIO
import matplotlib.pyplot as plt
import h5py


#f = h5py.File('convergence_CENTRED.h5', 'r')

#print(f.keys())
#print(f['outputParams'].keys())


#####################################
# 3. Do stuff with the finished scan
#####################################

csp_CENTRED = DREAMIO.LoadHDF5AsDict('convergence_CENTRED.h5')
csp_QUICK   = DREAMIO.LoadHDF5AsDict('convergence_QUICK.h5')
csp_UPWIND  = DREAMIO.LoadHDF5AsDict('convergence_UPWIND.h5')
csp_UPWIND2 = DREAMIO.LoadHDF5AsDict('convergence_UPWIND2.h5')
csp_SMART   = DREAMIO.LoadHDF5AsDict('convergence_SMART.h5')
csp_MUSCL   = DREAMIO.LoadHDF5AsDict('convergence_MUSCL.h5')
csp_OSPRE   = DREAMIO.LoadHDF5AsDict('convergence_OSPRE.h5')
csp_TCDF   = DREAMIO.LoadHDF5AsDict('convergence_TCDF.h5')

nps  = csp_CENTRED['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['scanval']
nxis = csp_CENTRED['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['scanval']

print(nps)
print(nxis)

rr_CENTRED_pscan   = csp_CENTRED['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_CENTRED_xiscan  = csp_CENTRED['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_QUICK_pscan     = csp_QUICK['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_QUICK_xiscan    = csp_QUICK['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_UPWIND_pscan    = csp_UPWIND['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_UPWIND_xiscan   = csp_UPWIND['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_UPWIND2_pscan   = csp_UPWIND2['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_UPWIND2_xiscan  = csp_UPWIND2['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_SMART_pscan     = csp_SMART['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_SMART_xiscan    = csp_SMART['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_MUSCL_pscan     = csp_MUSCL['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_MUSCL_xiscan    = csp_MUSCL['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_OSPRE_pscan     = csp_OSPRE['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_OSPRE_xiscan    = csp_OSPRE['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_TCDF_pscan      = csp_TCDF['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_TCDF_xiscan     = csp_TCDF['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']


yMaxPlot = 7e19
fig, axes = plt.subplots(1,2)

axes[0].plot(nps, rr_CENTRED_pscan,  label='Centered', linewidth=2, color='black', linestyle='dashed')
axes[0].plot(nps, rr_QUICK_pscan,    label='QUICK',    linewidth=2, color='blue',  linestyle='dashed')
axes[0].plot(nps, rr_UPWIND_pscan,   label='UPWIND',   linewidth=2, color='red',   linestyle='dashed')
axes[0].plot(nps, rr_UPWIND2_pscan,  label='UPWIND-2', linewidth=2, color='green', linestyle='dashed')
axes[0].plot(nps, rr_SMART_pscan,    label='SMART',    linewidth=2, color='black', linestyle='solid')
axes[0].plot(nps, rr_MUSCL_pscan,    label='MUSCL',    linewidth=2, color='blue',  linestyle='solid')
axes[0].plot(nps, rr_OSPRE_pscan,    label='OSPRE',    linewidth=2, color='red',   linestyle='solid')
axes[0].plot(nps, rr_TCDF_pscan,     label='TCDF',     linewidth=2, color='green', linestyle='solid')
axes[0].set_xlabel("$N_\\mathrm{p}$")
axes[0].set_ylabel("Runaway rate")
axes[0].set_xlim(left=nps[0], right=nps[-1])
axes[0].set_ylim(bottom=0, top=yMaxPlot)


axes[1].plot(nxis, rr_CENTRED_xiscan,  label='Centered', linewidth=2, color='black',  linestyle='dashed')
axes[1].plot(nxis, rr_QUICK_xiscan,    label='QUICK',    linewidth=2, color='blue',   linestyle='dashed')
axes[1].plot(nxis, rr_UPWIND_xiscan,   label='UPWIND',   linewidth=2, color='red',    linestyle='dashed')
axes[1].plot(nxis, rr_UPWIND2_xiscan,  label='UPWIND-2', linewidth=2, color='green',  linestyle='dashed')
axes[1].plot(nxis, rr_SMART_xiscan,    label='SMART',    linewidth=2, color='black',  linestyle='solid')
axes[1].plot(nxis, rr_MUSCL_xiscan,    label='MUSCL',    linewidth=2, color='blue',   linestyle='solid')
axes[1].plot(nxis, rr_OSPRE_xiscan,    label='OSPRE',    linewidth=2, color='red',    linestyle='solid')
axes[1].plot(nxis, rr_TCDF_xiscan,     label='TCDF',     linewidth=2, color='green',    linestyle='solid')
axes[1].set_xlabel("$N_\\mathrm{xi}$")
axes[1].set_ylabel("Runaway rate")
axes[1].set_xlim(left=nxis[0], right=nxis[-1])
axes[1].set_ylim(bottom=0, top=yMaxPlot)

axes[1].legend()
plt.show()

pMax = 0.3
pTe = np.sqrt(2*100/511e3)
dp = pMax/(pTe*nps)
dxi = 2/nxis

