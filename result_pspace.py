# -*- coding: utf-8 -*-
"""
Surface plot of key rate for scissors parameter space

Created on Tue May  7 15:02:55 2019

@author: Eduardo Villasenor 
"""

import src.cv_system as cv
import src.measurements as measurements
import numpy as np

#from mayavi import mlab

############################################ CALCULATIONS

## Parameters
N = 3
mpne = 0.001
f = 0.95
option = 'none'
eta = 0.01


theta = np.arccos(np.sqrt(eta))
r_eve = np.arcsinh(np.sqrt(mpne))

## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)

key_rates = []
ps = []
ks = np.linspace(0000.1, .999, 20)
mus = np.linspace(0000.1, 1.5, 20)

for k in ks:
    k_temp = []
    p_temp = []
    for mu in mus:
        print("--->", k, mu)
        sys.reset_state(2)
        r = np.arcsinh(np.sqrt(mu))
        sys.apply_TMS(r, [0, 1])
                
    
        # Transmitter Scissors
        if option == 'tsc':
            p_success = sys.apply_scissors_exact(k, 1)
            print("P SUCCESS:", p_success)
        
        if option == 'tps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)
        
        if option == 'none':
            p_success = 1
    
        sys.add_TMSV(r_eve)

    
        sys.apply_BS(theta, [1, 2])
    
    
        # Receiver Scissors
        if option == 'rsc':
            p_success = sys.apply_scissors_exact(k, 1)
            print("P SUCCESS:", p_success) 
            
            
        if option == 'rps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)

        kr = measurements.key_rate(sys, f, p_success)
        print("Key rate:", kr)
        k_temp += [kr]
        p_temp += [p_success]
    key_rates += [k_temp]
    ps += [p_temp]


# Save the resuls
filename = "data/result_SCpspace_" + option 
key = np.array(key_rates)
np.save(filename, key_rates)

filename = "data/result_SCpspace_p_" + option 
key_rates = np.array(ps)
np.save(filename, ps)

filename_ind1 = "data/indeces_SCpspace_k_" + option 
np.save(filename_ind1, ks)

filename_ind2 = "data/indeces_SCpspace_m_" + option 
np.save(filename_ind2, mus)


############################################ PLOT

import plot_pspace as plt

plt.plot(option)

