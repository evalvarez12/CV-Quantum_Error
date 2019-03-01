# -*- coding: utf-8 -*-
"""
Operations to realize entanglement distillation via quantum scisors

Created on Thu Feb 28 11:04:31 2019

@author: Eduardo Villasenor
"""

import numpy as np
import qutip as qt
import operations as ops
import tools

def scissor(rho_in, kappa):
    N = rho_in.dims[0][0]
    b = qt.destroy(N)
    
    theta1 = np.arccos(np.sqrt(kappa))
    theta2 = np.pi/4

    # Operator space order:
    # c -> c id id
    # b -> id b id
    # a -> id id a 
    c, b, a = ops.tritter([b, b, b], theta1, theta2)
    
    proyector0 = qt.basis(N)
    proyector1 = qt.basis(N, 1)
    proyector10 = qt.tensor([proyector0, qt.identity(N), proyector1])
    
    vacuum = qt.basis(N) * qt.basis(N).dag()
    vacuum_b = tools.tensor(vacuum, [N]*2, 1)
    vacuum_c = tools.tensor(vacuum, [N]*2, 0)
        
    
    print(proyector10.dims)
    print(vacuum_b.dims)
    print(vacuum_c.dims)
    print(c.dims)
    a_out = proyector10.dag()* rho_in * vacuum_b * c.dag() * vacuum_c * c * proyector10
    return a_out


N =5 
kappa = 0.5

a = qt.destroy(N)

print(scissor(a, kappa))