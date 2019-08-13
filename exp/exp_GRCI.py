# -*- coding: utf-8 -*-
"""
Testing GRCI results

Created on Wed Mar  6 16:19:53 2019

@author: Eduardo Villasenor
"""

import cv_system as cv
import numpy as np
import measurements
import matplotlib.pyplot as plt

# Initial parameters
N = 20

# Sweep parameters
kappa = 0.005
mu = 1.8
r = np.arcsinh(np.sqrt(mu))
print("r:", r)
# Fixed parameters
#mu_aux = 0.01
#r_aux = np.arcsinh(np.sqrt(mu_aux))


sys = cv.System(N, Nmodes=2, cm=False)
sys.apply_TMS(r, pos=[0, 1])
sys.save_state()

results = []
etas = np.linspace(0.01, .99, 20)
for eta in etas:
    sys.apply_loss_channel(eta, 1)
    rci = measurements.RCI(sys, [0])
    sys.load_state()
    print(eta, rci)
    results += [rci]


plt.plot(etas, results, 'bo')
ref = -np.log2(1 - etas)
plt.plot(etas, ref, 'ko')
plt.show()
