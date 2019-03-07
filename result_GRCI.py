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
N = 5

# Sweep parameters
kappa = 0.005
mu = 0.05

# Fixed parameters
mu_aux = 0.01
eta = 0.01

vacuum = qt.basis(N) * qt.basis(N).dag()
vacuum = qt.tensor(vacuum, vacuum)

S =  ops.tmsqueeze(N, mu)

rho_in = S * vacuum * S.dag()

rho = bs.loss_channel_applyU(rho_in, 0, eta)

rho = scissors.scissor_1NLA(rho, kappa, mu_aux)


rci = ops.RCI(rho, 1)
print(rci)