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
N = 5

# Fixed parameters
mu_aux = 0.1
r_aux = np.arcsinh(np.sqrt(mu_aux))
eta = 0.01


sys1 = cv.System(N, Nmodes=2, cm=False)

sys2 = cv.System(N, Nmodes=2, cm=False)

sys1.save_state()
results_sc = []
results_ps = []
probabilities_sc = []
probabilities_ps = []

for kappa in np.linspace(0.0, 0.01, 20):
#for kappa in np.linspace(0.0001, 0.03, 20):
    res_sc = []
    res_ps = []
    ps_sc = []
    ps_ps = []
    for mu in np.linspace(0.0, .1, 20):
#    for mu in np.linspace(0.0001, .6, 20):

        r = np.arcsinh(np.sqrt(mu))
        sys1.apply_TMS(r, pos=[0, 1])


        sys1.apply_loss_channel(eta, 1)

#        sys2.set_state(sys1.state)


        # Quantum scissors
        # p_sc = sys1.apply_scissors(kappa, r_aux, 1)
        p_sc = sys1.apply_scissors_options(kappa, r_aux, 1, 'b')
#        p_sc = sys1.apply_scissors_exact(kappa, 1)

#        p_sc = sys1.apply_scissors_options(kappa, r_aux, 1, 'b')

        # Photon substraction
#        p_ps = sys2.apply_photon_subtraction(1- kappa, 1)
        p_ps = 1

#        print("Reversed")
        rci_sc = measurements.CI(sys1, [0])
#        print("Coherent")
#        ci = measurements.CI(sys2, [1])

#        rci_ps = measurements.CI(sys2, [0])
        rci_ps = 0

        print(mu, kappa, p_sc, rci_sc)
#
#        print(rci > ci)
#        if rci > ci:
#            marker += 1

        sys1.load_state()

        res_sc += [rci_sc]
        res_ps += [rci_ps]
        ps_sc += [p_sc]
        ps_ps += [p_ps]

    results_sc += [res_sc]
    results_ps += [res_ps]
    probabilities_sc += [ps_sc]
    probabilities_ps += [ps_ps]

filename_sc = "data/rci_plot_SC"

results_sc = np.array(results_sc)
np.save(filename_sc, results_sc)

filename_p_sc = "data/p_plot_SC"

probabilities_sc = np.array(probabilities_sc)
np.save(filename_p_sc, probabilities_sc)



filename_ps = "data/rci_plot_PS"

results_ps = np.array(results_ps)
np.save(filename_ps, results_ps)

filename_p_ps = "data/p_plot_PS"

probabilities_ps = np.array(probabilities_ps)
np.save(filename_p_ps, probabilities_ps)

#################

import plt_GRCI as plt

plt.plot()
