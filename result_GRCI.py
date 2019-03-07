# -*- coding: utf-8 -*-
"""
Replicating Gaussian RCI plots

Created on Wed Mar  6 16:19:53 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import operations as ops
import scissors
import beam_splitter as bs

# Initial parameters
N = 10

# Sweep parameters
kappa = 0.005
mu = 0.05

# Fixed parameters
mu_aux = 0.01
eta = 0.01

vacuum = qt.basis(N) * qt.basis(N).dag()
vacuum = qt.tensor(vacuum, vacuum)

results = []

for kappa in np.linspace(0.001, 0.01, 20):
    r2 = []
    for mu in np.linspace(0.001, .1, 20):

        S =  ops.tmsqueeze(N, mu)
        
        rho_in = S * vacuum * S.dag()
        
        rho = bs.loss_channel_applyU(rho_in, 0, eta)
        
        rho = scissors.scissor_1NLA(rho, kappa, mu_aux)
        
        
        rci = ops.RCI(rho, 1)

        r2 +=[rci]
    results += [r2]



filename = "rci_plot_1NLA"

results = np.array(results)
np.save(filename, results)
