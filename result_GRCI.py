# -*- coding: utf-8 -*-
"""
Replicating Gaussian RCI plots

Created on Wed Mar  6 16:19:53 2019

@author: Eduardo Villasenor
"""

import src.cv_system as cv
import src.measurements as measurements
import numpy as np


################################## CALCULATIONS

# Initial parameters
N = 20
mu_aux = 0.1
r_aux = np.arcsinh(np.sqrt(mu_aux))
eta = 0.01
option = 'rct'

sys = cv.System(N, Nmodes=2, cm=False)


sys.save_state()
rcis = []
ps = []
probabilities_sc = []
probabilities_ps = []

ks = np.linspace(0.0, 0.01, 20)
mus = np.linspace(0.0, .1, 20)

for k in ks:
#for kappa in np.linspace(0.0001, 0.03, 20):
    rci_temp = []
    p_temp = []
    for mu in mus:
#    for mu in np.linspace(0.0001, .6, 20):

        r = np.arcsinh(np.sqrt(mu))
        sys.apply_TMS(r, pos=[0, 1])


        sys.apply_loss_channel(eta, 1)


        # Quantum scissors
        if option == 'rsc':
            p = sys.apply_scissors_exact(k, 1)
#            p = sys1.apply_scissors(kappa, r_aux, 1)

        # Photon substraction
        elif option == 'rps':
           p = sys.apply_photon_subtraction(k, 1)
        
        elif option == 'rct':
           p = sys.apply_photon_catalysis(1, k, 1) 
        
        elif option == 'none':
           p = 1

        rci = measurements.CI(sys, [0])

        print(mu, k, p, rci)
        sys.load_state()

        rci_temp += [rci]
        p_temp += [p]

    rcis += [rci_temp]
    ps += [p_temp]

filename_rci = "data/result_GRCI_" + option
rcis = np.array(rcis)
np.save(filename_rci, rcis)

filename_p = "data/result_GRCI_p_" + option
ps = np.array(ps)
np.save(filename_p, ps)


filename_ind1 = "data/indeces_GRCI_k_" + option 
np.save(filename_ind1, ks)

filename_ind2 = "data/indeces_GRCI_m_" + option 
np.save(filename_ind2, mus)

#################

import plot_GRCI as plt

plt.plot(option)
