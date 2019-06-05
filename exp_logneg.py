# -*- coding: utf-8 -*-
"""
Surface plot of key rate for scissors parameter space

Created on Tue May  7 15:02:55 2019

@author: Eduardo Villasenor 
"""

import src.cv_system as cv
import src.measurements as measurements
import numpy as np
import matplotlib.pyplot as plt


############################################ CALCULATIONS

## Parameters
N = 30
k = 0.95
option = 'dual'

## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)

els = []
ps = []
mus = np.linspace(0.001, 2, 20)

for mu in mus:
    print("--->", mu)
    sys.reset_state(2)

#    r = np.arcsinh(np.sqrt(mu))
    r = mu
    sys.apply_TMS(r, [0, 1])
        
    if option == 'none':
        p_success = 1
    
    elif option == 'single':
        p_success = sys.apply_photon_subtraction(k, 1)
        print("P SUCCESS:", p_success)

    elif option == 'dual':
        p_success0 = sys.apply_photon_subtraction(k, 0)
        p_success1 = sys.apply_photon_subtraction(k, 1)
        p_success = p_success0 * p_success1
        print("P SUCCESS:", p_success)
            
    el = measurements.log_neg(sys.state, [1, 0])
#    el = measurements.log_neg_Gauss(sys, [1, 0])

    print("Logarithmic Negativity", el, p_success)
    els += [el]
    ps += [p_success]


# Save the resuls
filename = "data/exp_logneg_el_" + option 
els = np.array(els)
np.save(filename, els)

filename = "data/exp_logneg_p_" + option 
key_rates = np.array(ps)
np.save(filename, ps)

filename_ind1 = "data/indeces_logneg_mu_" + option 
np.save(filename_ind1, mus)



############################################ PLOT

list_plots = ['none', 'single', 'dual']

fig = plt.figure()
indeces = []
els = []

for ext in list_plots :

    filename = "data/exp_logneg_el_" + ext 
    el = np.load(filename + '.npy')
    
    filename = "data/exp_logneg_p_" + ext
    ps = np.load(filename + '.npy')
    
    filename_ind1 = "data/indeces_logneg_mu_" + ext 
    inds = np.load(filename_ind1 + '.npy')

    indeces += [inds]
    els += [el]


lines_types = ['k*-', 'b*-', 'r*-']
for i in range(len(list_plots)):
    plt.plot(indeces[i], els[i], lines_types[i])

# Plot the surface.

plt.title(r"Logarithmic Negativity")
plt.xlabel(r'$\mu$', size=15)
plt.ylabel(r'$E_N$', size=10)

#fig = plt.figure()
#plt.plot(inds, ps)
#plt.title(r"P success Scissors")
#plt.xlabel(r'$p_{success}$', size=15)
#plt.ylabel(r'Key rate', size=10)

plt.show()
