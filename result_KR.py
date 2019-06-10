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

#from mayavi import mlab

############################################ CALCULATIONS

## Parameters
N = 20
mpne = 0.001
f = 0.95
option = 'rct'
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
        
        elif option == 'tps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)
        
        elif option == 'none':
            p_success = 1
    
        elif option == 'tct':
            p_success = sys.apply_photon_catalysis(1, k, 1)
            print("P_SUCCESS:", p_success)
    
        sys.add_TMSV(r_eve)

    
        sys.apply_BS(theta, [1, 2])
    
    
        # Receiver Scissors
        if option == 'rsc':
            p_success = sys.apply_scissors_exact(k, 1)
            print("P SUCCESS:", p_success) 
            
        elif option == 'rps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)

        elif option == 'rct':
            p_success = sys.apply_photon_catalysis(1, k, 1)
            print("P_SUCCESS:", p_success)
            
        kr = measurements.key_rate(sys, f, p_success)
        print("Key rate:", kr)
        k_temp += [kr]
        p_temp += [p_success]
    key_rates += [k_temp]
    ps += [p_temp]



# File name parameters
k_name = 'var'
mu_name = 'var'
eta_name = eta
measurement = "KR"
measurementp = "KR_p"


# Save the resuls
filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
key = np.array(key_rates)
np.save(filename, key_rates)

filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
key_rates = np.array(ps)
np.save(filenamep, ps)

filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
np.save(filename_ind1, ks)

filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='mu')
np.save(filename_ind2, mus)


############################################ PLOT

import plots 

plots.KR(option, N, eta)

