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
N = 2
mpne = 0.001
mu = .01
f = 0.95
option = 'rct'

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
        
        elif option == 'tps':
            p_success = sys.apply_photon_subtraction(k, 1)
            print("P SUCCESS:", p_success)
        
        elif option == 'none':
            p_success = 1
    
        sys.add_TMSV(r_eve)

        theta = np.arccos(np.sqrt(eta))
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
            print("P SUCCESS:", p_success)

        kr = measurements.key_rate(sys, f, p_success)
        print("Key rate:", kr)
        k_temp += [kr]
        p_temp += [p_success]
    key_rates += [k_temp]
    ps += [p_temp]


# File name parameters
k_name = 'var'
mu_name = mu
eta_name = 'var'
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

filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='eta')
np.save(filename_ind2, etas)


############################################ PLOT

import plots as plt

plt.KR2(option, N, mu)