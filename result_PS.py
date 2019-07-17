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
options = ['none', 'tps', 'rps', 'tqs', 'rqs']
#options = ['none', 'tps', 'rps']

#options = ['none', 'tqs', 'rqs']

for option in options:
#for i in [0]:
#    option = 'rsc'

    # Parameters
    N = 20
    r = .92
    r_eve = 0.033

    print("Protocol:", option)
    if option == 'rqs':
        N = 10
    
    # Operations options
    k_ps = 0.95
    k_qs = 0.05
    
    
    
    ## Initialize state
    sys = cv.System(N, Nmodes=2, cm=False)
    sys.apply_TMS(r, [0, 1])
    # BAD TMSV
    #sys.replace_current_state_w_bad_TMSV(mean_photon_number)
    
    
    # Transmitter Photon subtraction
    if option == 'tps':
        p_success = sys.apply_photon_subtraction(k_ps, 1)
        print("P SUCCESS:", p_success)
    
    # No Photon subtraction
    elif option == 'none':
        p_success = 1
    
    # Transmitter Scissors Exact
    elif option == 'tqs':
        p_success = sys.apply_scissors_exact(k_qs, 1)
        print("P SUCCESS:", p_success)
    
    # Evesdropper collective attack
    sys.add_TMSV(r_eve)

    
    sys.set_quadratures_basis()
    
    
    # Save current state of the system
    sys.save_state()
    
    key_rates = []
    
    ##### For r=.92
    if option == 'none' or option == 'tps':
        tes = np.logspace(-2, 0, base=10, num=32)
    
    if option == 'rps':
        tes = np.logspace(-3, 0, base=10, num=38)
        tes = tes[6:]
        
    if option == 'rqs' or option == 'tqs':   
        tes = np.linspace(1, 0.98, num=3)
        
    
    
    ##### For r=.12
#    if option == 'none' or option == 'tps':
#        tes = np.logspace(-1, 0, base=10, num=35)[10:]
#    
#    if option == 'rps':
#        tes = np.logspace(-2, 0, base=10, num=35)[10:]
##        tes = tes[6:]
#    
#    if option == 'rqs':
##        tes = np.logspace(-1, 0, base=10, num=3)
#        tes = np.logspace(-2, 0, base=10, num=40)[12:]
#        
#    if option == 'tqs':
#        tes = np.logspace(-2, 0, base=10, num=40)[12:]
        
    #tes = np.linspace(.9, 1, 10)
    #tes = [1.]
    
    for te in tes:
        print("->", te)
        sys.load_state()
    
        theta = np.arccos(np.sqrt(te))
        sys.apply_BS(theta, [1, 2])
    
        # Receiver Photon subtraction
        if option == 'rps':
            p_success = sys.apply_photon_subtraction(k_ps, 1)
            print("P SUCCESS:", p_success)
            
        # Receiver Scissors Exact
        elif option == 'rqs':
            p_success = sys.apply_scissors_exact(k_qs, 1)
            sys.state = sys.state.permute([2,3,1,0])
            print("P SUCCESS:", p_success)
    
    
        kr = measurements.key_rate(sys, 1, p_success)
        key_rates += [kr]
        print("KR:", kr)
    
    
    params = ["r=" + str(r), "r_eve=" + str(r_eve), "k_ps=" + str(k_ps),
              "k_qs=" + str(k_qs)]
    # Save the resuls
    filename = names.measurements_line(N, 'KR2', params, option)
    key_rates = np.array(key_rates)
    #print(key_rates)
    np.save(filename, key_rates)
    print(f"-------------")
    print(filename)
    filename_ind = names.indeces_line(N, 'KR2', params, option, 'eta')
    np.save(filename_ind, tes)


############################################ PLOT
import plot_PS as plot

plot.plot2(N, params)