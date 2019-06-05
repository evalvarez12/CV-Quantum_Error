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
N = 20
option = 'rct'
eta = 0.01


## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)

els = []
ps = []
ks = np.linspace(0000.1, .999, 20)
mus = np.linspace(0000.1, 1.5, 20)

for k in ks:
    el_temp = []
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
        
        elif option == 'tps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)
        
        elif option == 'none':
            p_success = 1
    
        elif option == 'tct':
            p_success = sys.apply_photon_catalysis(1, k, 1)
            print("P SUCCESS:", p_success)
    
    
        sys.apply_loss_channel(eta, 1)
    
    
        # Receiver Scissors
        if option == 'rsc':
            p_success = sys.apply_scissors_exact(k, 1)
            print("P SUCCESS:", p_success) 
            
        elif option == 'rps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)
            
        elif option == 'rct':
            p_success = sys.apply_photon_catalysis(1, k, 1)
            print("P SUCCESS:", p_success)

        el = measurements.log_neg(sys.state, [0, 1])
        print("Logarithmic Negativity:", el)
        el_temp += [el]
        p_temp += [p_success]
    els += [el_temp]
    ps += [p_temp]


# Save the resuls
filename = "data/result_EN_" + option 
els = np.array(els)
np.save(filename, els)

filename = "data/result_EN_p_" + option 
key_rates = np.array(ps)
np.save(filename, ps)

filename_ind1 = "data/indeces_EN_k_" + option 
np.save(filename_ind1, ks)

filename_ind2 = "data/indeces_EN_m_" + option 
np.save(filename_ind2, mus)


############################################ PLOT

import plot_EN as plt

plt.plot(option)
