# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 11:42:52 2019

@author: Eduardo Villasenor 
"""

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
options = ['none', 'tps', 'rps', 'tsc', 'rsc']
#options = ['none', 'tps', 'rps']


N = 10
mpn = 0.1
print("Protocol:", option)

# Operations options
k_ps = 0.5
k_sc = 0.01
k_ct = 0.5

r = np.arcsinh(np.sqrt(mpn))

etas = np.logspace(-2, 0, base=10, num=20)


for option in options:
#for i in [0]:
    
    ## Initialize state
    sys = cv.System(N, Nmodes=2, cm=False)
    sys.apply_TMS(r, [0, 1])
    
    
    # Transmitter Photon subtraction
    if option == 'tps':
        p_success = sys.apply_photon_subtraction(k_ps, 1)
        print("P SUCCESS:", p_success)
    
    elif option == 'tct':
        p_success = sys.apply_photon_catalysis(1, k_ct, 1)
        print("P SUCCESS:", p_success)
    
    # No Photon subtraction
    elif option == 'none':
        p_success = 1
    
    # Transmitter Scissors 
    elif option == 'tsc':
        p_success = sys.apply_scissors_exact(k_sc, 1)
        print("P SUCCESS:", p_success)
    #    print(sys.state)

    sys.save_state()
    
    els = []
    ps = []
    #tes = np.linspace(.9, 1, 10)
    #tes = [1.]
    
    for eta in etas:
        sys.load_state()
    
        sys.apply_loss_channel(eta, 1)
    
    
        # Receiver Scissors
        if option == 'rsc':
            p_success = sys.apply_scissors_exact(k_sc, 1)
            print("P SUCCESS:", p_success) 
            
        elif option == 'rps':
            p_success = sys.apply_photon_subtraction(k_ps, 1)
            print("P SUCCESS:", p_success)
            
        elif option == 'rct':
            p_success = sys.apply_photon_catalysis(1, k_ct, 1)
            print("P SUCCESS:", p_success)

        el = measurements.log_neg(sys.state, [0, 1])
        print("Logarithmic Negativity:", el)
        els += [el]
        ps += [p_success]
        

    
    params = ["mpn=" + str(mpn), "k_ps=" + str(k_ps),
              "k_sc=" + str(k_sc), "k_ct=" + str(k_ct)]
    # Save the resuls
    filename = names.measurements_line(N, 'EN2', params, option)
    els = np.array(els)
    #print(key_rates)
    np.save(filename, els)
    
    
    filename = names.measurements_line(N, 'EN2_p', params, option)
    ps = np.array(ps)
    #print(key_rates)
    np.save(filename, ps)
    
    filename_ind = names.indeces_line(N, 'EN2', params, option, 'eta')
    np.save(filename_ind, etas)


############################################ PLOT
import plot_EN2 as plot

plot.plot(N, params)