# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:36:08 2019

@author: Eduardo Villasenor 
"""

import qutip as qt
import numpy as np
import operations as ops
import beam_splitter as bs


# Define some initail parameters
N = 3

#S = qt.squeeze(N, eps)
#D = qt.displace(N, 1)
vacuum = qt.basis(N)

state_zero = qt.tensor(vacuum, vacuum)

a = qt.destroy(N)
b = qt.destroy(N)

theta = np.pi/4
a,b = bs.beam_splitter([a, b], theta)

# Hong-Ou-Mandel interference
# Two photons enter both sides of a 50/50 bs
state = a.dag()*b.dag()*state_zero

print(state)

print(state * state.dag())


# Testing now using Schodringuer picture and unitary evolution
state_sh = qt.basis(N, 1) * qt.basis(N, 1).dag()
state_sh = qt.tensor(state_sh, state_sh)
print(state_sh)
state_sh = bs.beam_splitter_applyU(state_sh, theta)
print(state_sh)

