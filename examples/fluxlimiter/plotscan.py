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
csp_SMART   = DREAMIO.LoadHDF5AsDict('convergence_SMART.h5')
csp_MUSCL   = DREAMIO.LoadHDF5AsDict('convergence_MUSCL.h5')
#csp_SMART_PE = DREAMIO.LoadHDF5AsDict('convergence_SMART_PE.h5')
#csp_MUSCL_PE = DREAMIO.LoadHDF5AsDict('convergence_MUSCL_PE.h5')

nps  = csp_CENTRED['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['scanval']
nxis = csp_CENTRED['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['scanval']

rr_CENTRED_pscan   = csp_CENTRED['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_CENTRED_xiscan  = csp_CENTRED['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_QUICK_pscan     = csp_QUICK['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_QUICK_xiscan    = csp_QUICK['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_SMART_pscan     = csp_SMART['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_SMART_xiscan    = csp_SMART['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_MUSCL_pscan     = csp_MUSCL['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_MUSCL_xiscan    = csp_MUSCL['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
'''
rr_SMART_PE_pscan  = csp_SMART_PE['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_SMART_PE_xiscan = csp_SMART_PE['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
rr_MUSCL_PE_pscan  = csp_MUSCL_PE['result']['hottailgrid.pgrid.np']['other.fluid.runawayRate']['outval']
rr_MUSCL_PE_xiscan = csp_MUSCL_PE['result']['hottailgrid.xigrid.nxi']['other.fluid.runawayRate']['outval']
'''
fig, axes = plt.subplots(1,2)

yMaxPlot = 4e19

axes[0].plot(nps, rr_CENTRED_pscan,  label='Centered', linewidth=2, color='blue', linestyle='dashed')
axes[0].plot(nps, rr_QUICK_pscan,    label='QUICK', linewidth=2, color='red',linestyle='dashed')
axes[0].plot(nps, rr_SMART_pscan,    label='SMART', linewidth=2, color='green',linestyle='solid')
axes[0].plot(nps, rr_MUSCL_pscan,    label='MUSCL', linewidth=2, color='black',linestyle='solid')
#axes[0].plot(nps, rr_SMART_PE_pscan, label='SMART-PE', linewidth=2, color='green',linestyle='dotted')
#axes[0].plot(nps, rr_MUSCL_PE_pscan, label='MUSCL-PE', linewidth=2, color='black',linestyle='dotted')
axes[0].set_xlabel("$N_\\mathrm{p}$")
axes[0].set_ylabel("Runaway rate")
axes[0].set_xlim(left=nps[0], right=nps[-1])
axes[0].set_ylim(bottom=0, top=yMaxPlot)


axes[1].plot(nxis, rr_CENTRED_xiscan,  label='Centered',linewidth=2, color='blue', linestyle='dashed')
axes[1].plot(nxis, rr_QUICK_xiscan,    label='QUICK',linewidth=2,color='red',linestyle='dashed')
axes[1].plot(nxis, rr_SMART_xiscan,    label='SMART', linewidth=2, color='green',linestyle='solid')
axes[1].plot(nxis, rr_MUSCL_xiscan,    label='MUSCL', linewidth=2, color='black',linestyle='solid')
#axes[1].plot(nxis, rr_SMART_PE_xiscan, label='SMART-PE', linewidth=2, color='green',linestyle='dotted')
#axes[1].plot(nxis, rr_MUSCL_PE_xiscan, label='MUSCL-PE', linewidth=2, color='black',linestyle='dotted')
axes[1].set_xlabel("$N_\\mathrm{xi}$")
axes[1].set_ylabel("Runaway rate")
axes[1].set_xlim(left=nxis[0], right=nxis[-1])
axes[1].set_ylim(bottom=0, top=yMaxPlot)

axes[1].legend()
plt.show()



pMax = 2
pTe = np.sqrt(2*100/511e3)
dp = pMax/(pTe*nps)
dxi = 2/nxis

