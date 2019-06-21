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
N = 10
k_ps = 0.95
k_qs = 0.05
#option = 'single'

# eta = 0.315 = 5 dB
# eta = 0.1 = 10 dB
eta = 0.1
theta = np.arccos(np.sqrt(eta))

# Eve parameters
r_eve = 0.012

options = ['none', 'rps', 'tps', 'tqs', 'rqs']
## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)


mus = np.linspace(0.001, 1, 20)

for option in options:
    key_rates = []
    ps = []

    print("Options:", option)
    for mu in mus:
        print("--->", mu)
        sys.reset_state(2)
    
    #    r = np.arcsinh(np.sqrt(mu))
        r = mu
        sys.apply_TMS(r, [0, 1])
        
        
        # Transmitter Photon subtraction
        if option == 'tps':
            p_success = sys.apply_photon_subtraction(k_ps, 1)
            print("P SUCCESS:", p_success)
        
        # Transmitter Scissors 
        elif option == 'tqs':
            p_success = sys.apply_scissors_exact(k_qs, 1)
            print("P SUCCESS:", p_success)
        
        # No Photon subtraction
        elif option == 'none':
            p_success = 1
            
        sys.add_TMSV(r_eve)
            
        sys.apply_BS(theta, [1, 2])
        
        # Receiver Scissors
        if option == 'rqs':
            p_success = sys.apply_scissors_exact(k_qs, 1)
            print("P SUCCESS:", p_success) 
            
        elif option == 'rps':
            p_success = sys.apply_photon_subtraction(k_ps, 1)
            print("P SUCCESS:", p_success)


        sys.set_quadratures_basis()
        kr = measurements.key_rate(sys, 1, p_success)
        key_rates += [kr]
        print("KR:", kr)
        ps += [p_success]


    # Save the resuls
    filename = "data/res_KR_" + option 
    key_rates = np.array(key_rates)
    np.save(filename, key_rates)
    
    filename = "data/res_KR_p_" + option 
    key_rates = np.array(ps)
    np.save(filename, ps)
    
    filename_ind1 = "data/indeces_KR_r_" + option 
    np.save(filename_ind1, mus)



############################################ PLOT

list_plots = ['none', 'rps', 'tps', 'tqs', 'rqs']

fig = plt.figure()
indeces = []
els = []
pss = []

for ext in list_plots :

    filename = "data/res_KR_" + ext 
    el = np.load(filename + '.npy')
    
    filename = "data/res_KR_p_" + ext
    ps = np.load(filename + '.npy')
    
    filename_ind1 = "data/indeces_KR_r_" + ext 
    inds = np.load(filename_ind1 + '.npy')

    indeces += [inds]
    els += [el]
    pss += [ps]

lines_types = ['k*-', 'r*-', 'b*-', 'g*-', 'm*-']
legends = ['TMSV', 'rPS', 'tPS', 'tQS', 'rQS']
for i in range(len(list_plots)):
    x = -10*np.log10(np.exp(-2*indeces[i]))
#    y = np.multiply(ps[i], els[i])
    y = els[i]
    
    plt.plot(x, y, lines_types[i], label=legends[i])

# Plot the surface.

#plt.title(r"Logarithmic Negativity")
plt.title(r"Attenuation = 20 dB")
plt.xlabel(r'Squeezing (dB)', size=15)
plt.ylabel(r'$E_N$', size=15)

#plt.text(.05, .8, "Attenuation = 5 dB", size=15)
#fig = plt.figure()
#plt.plot(inds, ps)
#plt.title(r"P success Scissors")
#plt.xlabel(r'$p_{success}$', size=15)
#plt.ylabel(r'Key rate', size=10)

plt.legend(fontsize=15)

plt.show()
