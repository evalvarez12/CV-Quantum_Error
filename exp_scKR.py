# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 15:32:37 2019

@author: z5239621
"""


import src.cv_system as cv
import src.measurements as measurements
import numpy as np


############################################ CALCULATIONS
# Parameters
N = 10
mpn = 0.01
mpne = 0.001
f = 0.95
option = 'rsc'

print("N:", N)

te = 0.1

# Operations options
k_sc = 0.001

r = np.arcsinh(np.sqrt(mpn))
r_eve = np.arcsinh(np.sqrt(mpne))

## Initialize state
sys = cv.System(N, Nmodes=2, cm=False)
sys.apply_TMS(r, [0, 1])


# Evesdropper collective attack
sys.add_TMSV(r_eve)

sys.set_quadratures_basis()
tes = np.logspace(-2, 0, base=10, num=20)
#tes = np.linspace(.9, 1, 10)
#tes = [1.]


theta = np.arccos(np.sqrt(te))
sys.apply_BS(theta, [1, 2])


# Receiver Scissors Exact
p_success = sys.apply_scissors_exact(k_sc, 1)
print("P SUCCESS:", p_success)

kr = measurements.key_rate(sys, f, p_success)
print("--->", te, kr)

