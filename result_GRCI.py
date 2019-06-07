# -*- coding: utf-8 -*-
"""
Replicating Gaussian RCI plots

Created on Wed Mar  6 16:19:53 2019

@author: Eduardo Villasenor
"""

import src.cv_system as cv
import src.measurements as measurements
import src.names as names
import numpy as np


################################## CALCULATIONS

# Initial parameters
N = 2
mu_aux = 0.1
r_aux = np.arcsinh(np.sqrt(mu_aux))
eta = 0.01
option = 'none'

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



# File name parameters
k_name = 'var'
mu_name = 'var'
eta_name = eta
measurement = "GRCI"
measurementp = "GRCI_p"

# Save the resuls
filename = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option)
rcis = np.array(rcis)
np.save(filename, rcis)

filenamep = names.measurements(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurementp, protocol=option)
ps = np.array(ps)
np.save(filenamep, ps)

filename_ind1 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='k')
np.save(filename_ind1, ks)

filename_ind2 = names.indeces(N=N, eta=eta_name, k=k_name, mu=mu_name, measurement=measurement, protocol=option, index='mu')
np.save(filename_ind2, mus)

#################

import plots as plt

plt.GRCI(option, N, eta)
