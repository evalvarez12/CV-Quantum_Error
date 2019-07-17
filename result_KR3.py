# -*- coding: utf-8 -*-
"""
Replicaing Photon subtraction results of Minjgian's paper

Created on Mon Apr  1 14:12:46 2019

@author: Eduardo Villasenor
"""

import src.cv_system as cv
import src.measurements as measurements
import src.names as names
import numpy as np


############################################ CALCULATIONS
#options = ['none', 'rps', 'tqs']
options = ['rqs']
#options = ['none', 'tps', 'rps']

#options = ['none', 'tps', 'rps']


# Operations options
k_ps = 0.95
k_qs = 0.05

# Parameters
N = 10

r_eve = 0.033

for option in options:
    print(option)
#for i in [0]:
#    option = 'rsc'

    if option[1:] == 'ps':
        k = k_ps
        
    elif option[1:] == 'qs':
        k = k_qs

    else:
        k = 0
        
        
    print("Protocol:", option)
    if option == 'rqs':
        N = 10
    
    
    
    
    
    ## Initialize state
    sys = cv.System(N, Nmodes=2, cm=False)
    
    # BAD TMSV
    #sys.replace_current_state_w_bad_TMSV(mean_photon_number)
    
    
 
    
    
    # Save current state of the system
    sys.save_state()
    
    key_rates = []
    
    tes = np.logspace(-2, 0, base=10, num=50)
    rs = np.linspace(0.001, 0.93, num=50)
    
    for te in tes:
        for r in rs:
            print("->", te, r)
            
            sys.load_state()
            
            sys.apply_TMS(r, [0, 1])

            # Transmitter Photon subtraction
            if option == 'tps':
                p_success = sys.apply_photon_subtraction(k, 1)
                print("P SUCCESS:", p_success)
            
            # No Photon subtraction
            elif option == 'none':
                p_success = 1
            
            # Transmitter Scissors Exact
            elif option == 'tqs':
                p_success = sys.apply_scissors_exact(k, 1)
                print("P SUCCESS:", p_success)
            
            # Evesdropper collective attack
            sys.add_TMSV(r_eve)
            
            sys.set_quadratures_basis()
    
    
            theta = np.arccos(np.sqrt(te))
            sys.apply_BS(theta, [1, 2])
    
            # Receiver Photon subtraction
            if option == 'rps':
                p_success = sys.apply_photon_subtraction(k, 1)
                print("P SUCCESS:", p_success)
            
            # Receiver Scissors Exact
            elif option == 'rqs':
                p_success = sys.apply_scissors_exact(k, 1)
                print("P SUCCESS:", p_success)
    
    
            kr = measurements.key_rate(sys, 1, p_success)
            key_rates += [kr]
            print("KR:", kr)
    
    
    params = ["r=" + str(r), "r_eve=" + str(r_eve), "k=" + str(k)]
    # Save the resuls
    filename = names.measurements_line(N, 'KR3', params, option)
    key_rates = np.array(key_rates)
    #print(key_rates)
    np.save(filename, key_rates)
    
    filename_ind1 = names.indeces_line(N, 'KR3', params, option, 'eta')
    np.save(filename_ind1, tes)
    
    filename_ind2 = names.indeces_line(N, 'KR3', params, option, 'r')
    np.save(filename_ind2, rs)


############################################ PLOT
#import plot_PS as plot
#
#plot.plot2(N, params)