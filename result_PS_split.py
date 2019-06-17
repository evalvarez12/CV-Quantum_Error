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
options = ['none', 'tps', 'rps', 'rsc', 'tsc']
#options = ['none', 'tps', 'rps']

for option in options:
# for i in [0]:
    # Parameters
    # option = 'rsc'
    N = 10
    mpn = 0.01
    mpne = 0.005
    f = 0.95
    print("Protocol:", option)
    if option == 'rsc':
        N = 10

    # Operations options
    k_ps = 0.9
    k_sc = 0.1
    k_ct = 0.5

    # Modified Scissors options
    m_aux = .2

    r = np.arcsinh(np.sqrt(mpn))
    r_eve = np.arcsinh(np.sqrt(mpne))
    r_aux = np.arcsinh(np.sqrt(m_aux))

    ## Initialize state
    sys = cv.System(N, Nmodes=2, cm=False)
    sys.apply_TMS(r, [0, 1])
    # BAD TMSV
    #sys.replace_current_state_w_bad_TMSV(mean_photon_number)


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

    # Transmitter Scissors Exact
    elif option == 'tsc':
        p_success = sys.apply_scissors_exact(k_sc, 1)
        print("P SUCCESS:", p_success)
    #    print(sys.state)
    # Transmitter Scissors

    elif option == 'tsc_mod':
    #    p_success = sys.apply_scissors(k, r_aux, 1)
        p_success = sys.apply_scissors_options(k_sc, r_aux, 1, 'c')
        print("P SUCCESS:", p_success)
    #    print(sys.state)


    # Evesdropper collective attack
    sys.add_TMSV(r_eve)
    # BAD TMSV
    #sys.add_bad_TMSV(e_mpn)

    sys.set_quadratures_basis()


    # Save current state of the system
    sys.save_state()

    i_shareds = []
    i_stolens = []
    tes = np.logspace(-2, 0, base=10, num=20)
    #tes = np.linspace(.9, 1, 10)
    #tes = [1.]

    for te in tes:
        sys.load_state()

        theta = np.arccos(np.sqrt(te))
        sys.apply_BS(theta, [1, 2])

        # Receiver Photon subtraction
        if option == 'rps':
            p_success = sys.apply_photon_subtraction(k_ps, 1)
            print("P SUCCESS:", p_success)

        elif option == 'rct':
            p_success = sys.apply_photon_catalysis(1, k_ct, 1)
            print("P SUCCESS:", p_success)

        # Receiver Scissors Exact
        elif option == 'rsc':
            p_success = sys.apply_scissors_exact(k_sc, 1)
            print("P SUCCESS:", p_success)

        # Receiver Scissors
        elif option == 'rsc_mod':
            p_success = sys.apply_scissors(k_sc, r_aux, 1)
            print("P SUCCESS:", p_success)


    #    print(sys.cm)
    #    print(sys.get_full_CM())
        i_shared, i_stolen = measurements.key_rate_split(sys, f, p_success)
    #    kr = measurements.key_rate_compare(sys, f, p_success, mpn, mpne, te)
        i_shareds += [i_shared]
        i_stolens += [i_stolen]
        print("--->", te)


    params = ["mpn=" + str(mpn), "mpne=" + str(mpne), "f=" + str(f) , "k_ps=" + str(k_ps),
              "k_sc=" + str(k_sc), "k_ct=" + str(k_ct)]
    # Save the resuls
    filename1 = names.measurements_line(N, 'I_shared', params, option)
    i_shareds = np.array(i_shareds)
    np.save(filename1, i_shareds)

    filename2 = names.measurements_line(N, 'I_stolen', params, option)
    i_stolens = np.array(i_stolens)
    np.save(filename2, i_stolens)

    filename_ind = names.indeces_line(N, 'I_s', params, option, 'eta')
    np.save(filename_ind, tes)


############################################ PLOT
import plot_PS as plot

plot.plot_split(N, params)
