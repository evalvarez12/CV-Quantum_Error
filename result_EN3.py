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

### Parameters
#N = 20
#k_ps = 0.95
#k_qs = 0.05
##option = 'single'
#
## eta = 0.315 = 5 dB
## eta = 0.1 = 10 dB
#eta = .001
#
#options = ['none', 'rps', 'tps', 'rqs', 'tqs']
### Initialize system
#sys = cv.System(N, Nmodes=2, cm=False)
#
#
#mus = np.linspace(0.001, 1, 20)
#
#for option in options:
#    els = []
#    ps = []
#
#    print("Options:", option)
#    for mu in mus:
#        print("--->", mu)
#        sys.reset_state(2)
#    
#    #    r = np.arcsinh(np.sqrt(mu))
#        r = mu
#        sys.apply_TMS(r, [0, 1])
#        
#        
#        # Transmitter Photon subtraction
#        if option == 'tps':
#            p_success = sys.apply_photon_subtraction(k_ps, 1)
#            print("P SUCCESS:", p_success)
#        
#        # Transmitter Scissors 
#        elif option == 'tqs':
#            p_success = sys.apply_scissors_exact(k_qs, 1)
#            print("P SUCCESS:", p_success)
#        
#        # No Photon subtraction
#        elif option == 'none':
#            p_success = 1
#            
#        sys.apply_loss_channel(eta, 1)
#        
#        # Receiver Scissors
#        if option == 'rqs':
#            p_success = sys.apply_scissors_exact(k_qs, 1)
#            print("P SUCCESS:", p_success) 
#            
#        elif option == 'rps':
#            p_success = sys.apply_photon_subtraction(k_ps, 1)
#            print("P SUCCESS:", p_success)
#                
#        el = measurements.log_neg(sys.state, [1, 0])
#    #    el = measurements.log_neg_Gauss(sys, [1, 0])
#    
#        print("Logarithmic Negativity", el)
#        els += [el]
#        ps += [p_success]
#
#
#    # Save the resuls
#    filename = "data/res4_logneg_el_" + option 
#    els = np.array(els)
#    np.save(filename, els)
#    
#    filename = "data/res4_logneg_p_" + option 
#    key_rates = np.array(ps)
#    np.save(filename, ps)
#    
#    filename_ind1 = "data/indeces4_res_logneg_mu_" + option 
#    np.save(filename_ind1, mus)
##


############################################ PLOT

list_plots = ['none', 'tps', 'rps', 'tqs', 'rqs']

fig = plt.figure()
indeces = []
els = []
pss = []

for ext in list_plots :

    filename = "data/res_logneg_el_" + ext 
    el = np.load(filename + '.npy')
    
    filename = "data/res_logneg_p_" + ext
    ps = np.load(filename + '.npy')
    
    filename_ind1 = "data/indeces_res_logneg_mu_" + ext 
    inds = np.load(filename_ind1 + '.npy')

    indeces += [inds]
    els += [el]
    pss += [ps]

lines_types = ['k*-', 'b*-', 'r*-', 'g*-', 'm*-']
legends = ['TMSV', 'tPS', 'rPS', 'tQS', 'rQS']
for i in range(len(list_plots)):
    x = -10*np.log10(np.exp(-2*indeces[i]))
#    y = np.multiply(ps[i], els[i])
    y = els[i]
    
    plt.plot(x, y, lines_types[i], label=legends[i], linewidth=2)

# Plot the surface.

plt.xlim(-0.05, 8.8)
plt.ylim(-0.0, .9)

#plt.title(r"Logarithmic Negativity")
plt.rcParams["font.family"] = "Times New Roman"

plt.title(r"Attenuation = 5 dB", size=15)
plt.xlabel(r'Squeezing (dB)', size=15)
plt.ylabel(r'$E_N$', size=15)
plt.tick_params(axis='both', which='major', labelsize=15)

#plt.text(.05, .8, "Attenuation = 5 dB", size=15)
#fig = plt.figure()
#plt.plot(inds, ps)
#plt.title(r"P success Scissors")
#plt.xlabel(r'$p_{success}$', size=15)
#plt.ylabel(r'Key rate', size=10)

#plt.legend(fontsize=15)

plt.show()
