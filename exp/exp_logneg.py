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
sys.path.append("..") 


############################################ CALCULATIONS

## Parameters
N = 20
k = 0.95
#option = 'single'

eta = 0.01

options = ['none', 'single', 'dual']
## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)


mus = np.linspace(0.001, 1, 20)

for option in options:
    els = []
    ps = []

    print("Options:", option)
    for mu in mus:
        print("--->", mu)
        sys.reset_state(2)
    
    #    r = np.arcsinh(np.sqrt(mu))
        r = mu
        sys.apply_TMS(r, [0, 1])
            
#        sys.apply_loss_channel(eta, 1)
        
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
    
        print("Logarithmic Negativity", el)
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
pss = []

for ext in list_plots :

    filename = "data/exp_logneg_el_" + ext 
    el = np.load(filename + '.npy')
    
    filename = "data/exp_logneg_p_" + ext
    ps = np.load(filename + '.npy')
    
    filename_ind1 = "data/indeces_logneg_mu_" + ext 
    inds = np.load(filename_ind1 + '.npy')

    indeces += [inds]
    els += [el]
    pss += [ps]

lines_types = ['k*-', 'r*-', 'b*-']
for i in range(len(list_plots)):
    x = -10*np.log10(np.exp(-2*indeces[i]))
#    y = np.multiply(ps[i], els[i])
    y = els[i]
    
    plt.plot(x, y, lines_types[i])

# Plot the surface.

plt.title(r"Logarithmic Negativity")
plt.xlabel(r'Squeezing (dB)', size=15)
plt.ylabel(r'$E_N$', size=10)

#fig = plt.figure()
#plt.plot(inds, ps)
#plt.title(r"P success Scissors")
#plt.xlabel(r'$p_{success}$', size=15)
#plt.ylabel(r'Key rate', size=10)

plt.show()
