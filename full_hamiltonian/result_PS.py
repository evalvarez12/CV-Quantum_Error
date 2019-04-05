# -*- coding: utf-8 -*-
"""
Replicaing Photon subtraction results of Minjgian's paper

Created on Mon Apr  1 14:12:46 2019

@author: Eduardo Villasenor
"""

import cv_system as cv
import measurements
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.io


############################################ CALCULATIONS

## Parameters
#N = 20
#mean_photon_number = 1.3
#t = .9
#e_mpn = 0.001
#f = 0.95
#option = 'rps'
#
#
### Initialize state
#sys = cv.System(N, Nmodes=2)
#sys.apply_TMS(mean_photon_number, [0, 1])
## BAD TMSV
##sys.replace_current_state_w_bad_TMSV(mean_photon_number)
#
#
## Transmitter Photon subtraction
#if option == 'tps':
#    sys.add_vacuum()
#    sys.apply_BS(t, [1, 2])
#    p_success = sys.collapse_fock_state(1, 2)
#
## No Photon subtraction
#if option == 'nops':
#    p_success = 1
#
#
## Evesdropper collective attack
#sys.add_TMSV(e_mpn)
## BAD TMSV
##sys.add_bad_TMSV(e_mpn)
#
## Save current state of the system
#sys.save_state()
#
#key_rates = []
#tes =np.logspace(-10, 0, base=10)
#
#for te in tes:
#    sys.load_state()
#
#    sys.apply_BS(te, [1, 2])
#    
#    # Receiver Photon subtraction
#    if option == 'rps':
#        sys.add_vacuum()
#        sys.apply_BS(t, [1, 4])
#        p_success = sys.collapse_fock_state(1, 4)
#
#    key_rates += [measurements.key_rate(sys, f=f, p=p_success)]
#
#
## Save the resuls
#filename = "data/result_PS_" + option 
#key_rates = np.array(key_rates)
#
#np.save(filename, key_rates)
#


############################################ PLOT
filename1 = "data/result_PS_nops"
key_rates1 = np.load(filename1 + ".npy")
filename2 = "data/result_PS_tps"
key_rates2 = np.load(filename2 + ".npy")
filename3 = "data/result_PS_rps"
key_rates3 = np.load(filename3 + ".npy")

tes =np.logspace(-10, 0, base=10)

fig1, ax = plt.subplots()
ax.plot(tes, key_rates1, 'k-')
ax.plot(tes, key_rates2, 'b-')
ax.plot(tes, key_rates3, 'r-')


mg_ratesNO = scipy.io.loadmat('data/rate_no.mat')
mg_ratesR = scipy.io.loadmat('data/rate_r.mat')
mg_ratesT = scipy.io.loadmat('data/rate_t.mat')
mg_indNO = scipy.io.loadmat('data/ind_no.mat')
mg_indR = scipy.io.loadmat('data/ind_r.mat')
mg_indT = scipy.io.loadmat('data/ind_t.mat')

mg_ratesNO = mg_ratesNO['RateNo'][0]
mg_ratesR = mg_ratesR['RateR'][0]
mg_ratesT = mg_ratesT['RateT'][0]
mg_indNO = mg_indNO['indd'][0]
mg_indR = mg_indR['indR'][0]
mg_indT = mg_indT['indT'][0]

#mg_ind[mg_ind == 1] = 0
ax.plot(mg_indNO, mg_ratesNO, 'ko--')
ax.plot(mg_indR, mg_ratesR, 'ro--')
ax.plot(mg_indT, mg_ratesT, 'bo--')

#key_rates1b = np.load(filename1 + "BAD.npy")
#key_rates2b = np.load(filename2 + "BAD.npy")
#key_rates3b = np.load(filename3 + "BAD.npy")
#ax.plot(tes, key_rates1b, 'ko')
#ax.plot(tes, key_rates2b, 'bo')
#ax.plot(tes, key_rates3b, 'ro')

ax.set_yscale('log')


locmaj = matplotlib.ticker.LogLocator(base=10,numticks=10) 
ax.yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

#ax.minorticks_on()
plt.show()