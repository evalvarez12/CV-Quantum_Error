# -*- coding: utf-8 -*-
"""
Replicaing Photon subtraction results of Minjgian's paper

Created on Mon Apr  1 14:12:46 2019

@author: Eduardo Villasenor
"""

import cv_system as cv


### Parameters
N = 10
mean_photon_number = 1.3
t = 0.9
e_mpn = 0.001


# Initialize state
sys = cv.System(N, Nmodes=2)

sys.apply_TMS(mean_photon_number, [0, 1])


# Photon subtraction
sys.add_vacuum()
sys.apply_BS(t, [1, 2])
p_success = sys.collapse_fock_state(0, 2)

# Evesdropper collective attack
#sys.add_vacuum(2)
#sys.apply_TMS(e_mpn, [2, 3])
sys.add_TMSV(e_mpn)



cm = sys.get_full_CM()