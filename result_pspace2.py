# -*- coding: utf-8 -*-
"""
Surface plot of key rate for scissors parameter space

Created on Tue May  7 15:02:55 2019

@author: Eduardo Villasenor 
"""

import src.cv_system as cv
import src.measurements as measurements
import numpy as np

############################################ CALCULATIONS

## Parameters
N = 4
mpne = 0.001
mu = 1.2 
f = 0.95
option = 'none'

r = np.arcsinh(np.sqrt(mu))
r_eve = np.arcsinh(np.sqrt(mpne))

## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)
sys.apply_TMS(r, [0, 1])
sys.save_state()

key_rates = []
ps = []
etas = np.logspace(-3, -1, base=10, num=20)
ks = np.linspace(0000.1, .999, 20)

for k in ks:
    k_temp = []
    p_temp = []
    for eta in etas:
        print("--->", k, mu)
        sys.load_state()
        
    
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

        theta = np.arccos(np.sqrt(eta))
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
filename = "data/result_pspace2_" + option 
key_rates = np.array(key_rates)
np.save(filename, key_rates)

filename = "data/result_pspace2_p_" + option 
key_rates = np.array(ps)
np.save(filename, ps)

filename_ind1 = "data/indeces_pspace2_k_" + option 
np.save(filename_ind1, ks)

filename_ind2 = "data/indeces_pspace2_eta_" + option 
np.save(filename_ind2, etas)


############################################ PLOT

import plot_pspace2 as plt

plt.plot(option)