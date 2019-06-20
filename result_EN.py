# -*- coding: utf-8 -*-
"""
Surface plot of key rate for scissors parameter space

Created on Tue May  7 15:02:55 2019

@author: Eduardo Villasenor 
"""

import src.cv_system as cv
import src.measurements as measurements
import src.names as names
import numpy as np

############################################ CALCULATIONS

## Parameters
N = 10
option = 'rsc'
eta = 0.1


## Initialize system
sys = cv.System(N, Nmodes=2, cm=False)

els = []
ps = []
ks = np.linspace(0000.1, .999, 20)
mus = np.linspace(0.001, 1, 20)

for k in ks:
    el_temp = []
    p_temp = []
    for mu in mus:
        print("--->", k, mu)
        sys.reset_state(2)
        r = mu
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


# File name parameters
k_name = 'var'
mu_name = 'var'
eta_name = eta
measurement = "EN"
measurementp = "EN_p"

# Save the resuls
filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
els = np.array(els)
np.save(filename, els)

filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
ps = np.array(ps)
np.save(filenamep, ps)

filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
np.save(filename_ind1, ks)

filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='mu')
np.save(filename_ind2, mus)



############################################ PLOT

import plots as plt

plt.EN(option, N, eta)
