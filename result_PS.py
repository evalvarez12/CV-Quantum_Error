# -*- coding: utf-8 -*-
"""
Replicaing Photon subtraction results of Minjgian's paper

Created on Mon Apr  1 14:12:46 2019

@author: Eduardo Villasenor
"""

import cv_system as cv
import measurements
import numpy as np


############################################ CALCULATIONS

# Parameters
N = 20
mpn = 1.3
mpne = 0.001
f = 0.95
option = 'rps'

# Photon subtraction options
t = .9

# Scissors options
# Best  k=0.5 m_aux=0.6
k = .4
m_aux = .2



ps_theta = np.arccos(np.sqrt(t))
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
    p_success = sys.apply_photon_subtraction(t, 1)
    print("P SUCCESS:", p_success)

elif option == 'tct':
    p_success = sys.apply_photon_catalysis(1, t, 1)
    print("P SUCCESS:", p_success)

# No Photon subtraction
elif option == 'none':
    p_success = 1

# Transmitter Scissors
elif option == 'tsc':
#    p_success = sys.apply_scissors(k, r_aux, 1)
    p_success = sys.apply_scissors_options(k, r_aux, 1, 'c')
    print("P SUCCESS:", p_success)
#    print(sys.state)

# Transmitter Scissors Exact
elif option == 'tsc_e':
    p_success = sys.apply_scissors_exact(k, 1)
    print("P SUCCESS:", p_success)
#    print(sys.state)


# Evesdropper collective attack
sys.add_TMSV(r_eve)
# BAD TMSV
#sys.add_bad_TMSV(e_mpn)

sys.set_quadratures_basis()


# Save current state of the system
sys.save_state()

key_rates = []
tes = np.logspace(-2, 0, base=10, num=20)
#tes = np.linspace(.9, 1, 10)
#tes = [1.]

for te in tes:
    sys.load_state()

    theta = np.arccos(np.sqrt(te))
    sys.apply_BS(theta, [1, 2])

    # Receiver Photon subtraction
    if option == 'rps':
        p_success = sys.apply_photon_subtraction(t, 1)
        print("P SUCCESS:", p_success)
        
    elif option == 'rct':
        p_success = sys.apply_photon_catalysis(1, t, 1)
        print("P SUCCESS:", p_success)

    # Receiver Scissors
    elif option == 'rsc':
        p_success = sys.apply_scissors(k, r_aux, 1)
        print("P SUCCESS:", p_success)

    # Receiver Scissors Exact
    elif option == 'rsc_e':
        p_success = sys.apply_scissors_exact(k, 1)
        print("P SUCCESS:", p_success)

#    print(sys.cm)
#    print(sys.get_full_CM())
    kr = measurements.key_rate(sys, f, p_success)
#    kr = measurements.key_rate_compare(sys, f, p_success, mpn, mpne, te)
    key_rates += [kr]
    print("--->", te, kr)


# Save the resuls
filename = "data/result_PS_" + option
key_rates = np.array(key_rates)
#print(key_rates)
np.save(filename, key_rates)

filename_ind = "data/indeces_PS_" + option
np.save(filename_ind, tes)


############################################ PLOT
import plot_PS as plt

plt.plot()