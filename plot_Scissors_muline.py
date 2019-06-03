# -*- coding: utf-8 -*-
"""
Surface plot of key rate for scissors parameter space

Created on Tue May  7 15:02:55 2019

@author: Eduardo Villasenor 
"""

import cv_system as cv
import measurements
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


############################################ CALCULATIONS

## Parameters
N = 15
mpne = 0.001
f = 0.95
option = 'none'
eta = 0.5

k = 0.6

theta = np.arccos(np.sqrt(eta))
r_eve = np.arcsinh(np.sqrt(mpne))

## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)

key_rates = []
ps = []
mus = np.linspace(0.001, 2, 20)

for mu in mus:
    print("--->", mu)
    sys.reset_state(2)

    r = np.arcsinh(np.sqrt(mu))
    sys.apply_TMS(r, [0, 1])
        
    
    # Transmitter Scissors
    if option == 'tsc':
        p_success = sys.apply_scissors_exact(k, 1)
        print("P SUCCESS:", p_success)
    #    print(sys.state)
    
    if option == 'none':
        p_success = 1
    
    sys.add_TMSV(r_eve)

    
    sys.apply_BS(theta, [1, 2])
    
    
    # Receiver Scissors
    if option == 'rsc':
        p_success = sys.apply_scissors_exact(k, 1)
        print("P SUCCESS:", p_success) 
            

    kr = measurements.key_rate(sys, f, p_success)
    print("Key rate:", kr, p_success)
    key_rates += [kr]
    ps += [p_success]


# Save the resuls
filename = "data/result_SCe_pline_kr_" + option 
key_rates = np.array(key_rates)
np.save(filename, key_rates)

filename = "data/result_SCe_plinemu_p_" + option 
key_rates = np.array(ps)
np.save(filename, ps)

filename_ind1 = "data/indeces_SCe_plinemu_mu_" + option 
np.save(filename_ind1, mus)



############################################ PLOT


#option = 'tsc'

filename = "data/result_SCe_pline_kr_" + option 
krs = np.load(filename + '.npy')

filename = "data/result_SCe_plinemu_p_" + option
ps = np.load(filename + '.npy')

filename_ind1 = "data/indeces_SCe_plinemu_mu_" + option 
inds = np.load(filename_ind1 + '.npy')


# Plot the surface.

fig = plt.figure()
plt.plot(inds, krs)
plt.title(r"Transmitter side Scissors")
plt.xlabel(r'$\mu$', size=15)
plt.ylabel(r'Key rate', size=10)

#fig = plt.figure()
#plt.plot(inds, ps)
#plt.title(r"P success Scissors")
#plt.xlabel(r'$p_{success}$', size=15)
#plt.ylabel(r'Key rate', size=10)

plt.show()
