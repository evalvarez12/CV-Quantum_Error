# -*- coding: utf-8 -*-
"""
Fucntions recreate Braunstein CV-error correction

Created on Thu Feb 14 14:05:12 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import operations as ops
import tools 

# Define system paramters
N = 2
vacuum = qt.basis(N)
vacuum = qt.tensor([vacuum]*3)
# Define initial operators
a = qt.destroy(N)

# Encode
theta1 = np.arccos(1/np.sqrt(3))
theta2 = np.pi/4

a, b, c = ops.tritter([a]*3, theta1, theta2)


# Decode
a_out, b_out, c_out = ops.tritter([a, b, c], theta2, theta1, tensor=False)


# Homodyne measurement on ancillas



state_out = a_out.dag()*vacuum
print(state_out)


# Initial state for comparion
D = tools.tensor(qt.create(N), [N, N], 0)
state_compare = D*vacuum
print(state_compare)
