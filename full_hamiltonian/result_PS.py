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



############################################ CALCULATIONS

#### Parameters
N = 20
mean_photon_number = 1.3
t = 0.9
e_mpn = 0.001
f = 0.95


# Initialize state
sys = cv.System(N, Nmodes=2)
sys.apply_TMS(mean_photon_number, [0, 1])


# Photon subtraction
sys.add_vacuum()
sys.apply_BS(t, [1, 2])
p_success = sys.collapse_fock_state(0, 2)
print("p success:", p_success)

# Evesdropper collective attack
sys.add_TMSV(e_mpn)

# Save current state of the system
sys.save_state()

key_rates = []
tes =np.logspace(-10, 0, base=10)

for te in tes:
    sys.load_state()

    sys.apply_BS(te, [1, 2])
    key_rates += [measurements.key_rate(sys, f=f, p=1)]


# Save the resuls
filename = "data/result_PS_tps"
key_rates = np.array(key_rates)

np.save(filename, key_rates)



############################################ PLOT
filename1 = "data/result_PS_nops"
key_rates1 = np.load(filename + ".npy")
filename2 = "data/result_PS_tps"
key_rates2 = np.load(filename + ".npy")


tes =np.logspace(-10, 0, base=10)

fig1, ax = plt.subplots()
ax.plot(tes, key_rates1, 'b-')
ax.plot(tes, key_rates2, 'ro')
ax.set_yscale('log')


locmaj = matplotlib.ticker.LogLocator(base=10,numticks=10) 
ax.yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

#ax.minorticks_on()
plt.show()