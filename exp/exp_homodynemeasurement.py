# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 14:42:03 2019

@author: z5239621
"""

import qutip as qt
import numpy as np
import operations as ops


# Define some initail parameters
N = 10
eps = 1 + 1j


#D = qt.displace(N, 1)
S = qt.squeeze(N, eps)
vacuum = qt.basis(N)


state = S*vacuum

print("State", state)
print("<S>", ops.mean_homodyne(state, np.pi/2))


print("<S^2>", ops.var_homodyne(state, 0))
print("<S^2>", ops.var_homodyne(state, np.pi/2))

print(np.exp(-2*np.abs(eps))*N)