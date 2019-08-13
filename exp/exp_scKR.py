# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 15:32:37 2019

@author: z5239621
"""


import src.cv_system as cv
import src.measurements as measurements
import numpy as np
import qutip as qt

############################################ CALCULATIONS
# Parameters
N = 12
mu = 2
mpne = 2
f = 0.95
option = 'rsc'

te = 0.1

# Operations options
k_sc = 0.001

print("N:", N)
print("mu:", mu)
print("k:", k_sc)

r = np.arcsinh(np.sqrt(mu))
r_eve = np.arcsinh(np.sqrt(mpne))

## Initialize state
sys = cv.System(N, Nmodes=2, cm=False)
sys.apply_TMS(r, [0, 1])


# Evesdropper collective attack
sys.add_TMSV(r_eve)

state0 = sys.state

sys.set_quadratures_basis()
tes = np.logspace(-2, 0, base=10, num=20)
#tes = np.linspace(.9, 1, 10)
#tes = [1.]


theta = np.arccos(np.sqrt(te))
sys.apply_BS(theta, [1, 2])


state1 = sys.state

# Receiver Scissors Exact
p_success = sys.apply_scissors_exact(k_sc, 1)
print("P SUCCESS:", p_success)

kr = measurements.key_rate(sys, f, p_success)
print("--->", te, kr)

#print(state1)
#
#
#state0 = state0 * state0.dag()
#qt.hinton(state0, xlabels=[], ylabels=[])
#
#state1 = state1 * state1.dag()
#qt.hinton(state0, xlabels=[], ylabels=[])
##qt.hinton(state0)
#
#state2 = sys.state * sys.state.dag()
#qt.hinton(state2, xlabels=[], ylabels=[])
##qt.hinton(state)

