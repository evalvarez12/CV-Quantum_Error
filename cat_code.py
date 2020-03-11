# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 10:38:57 2020

@author: z5239621
"""

import qutip as qt
import numpy as np
import src.operations as ops
import src.wigner_plots as wp

#N = 30
#a = 3
#
#Np = 2 + 2 * np.exp(- 2* a**2)
#Nm = 2 - 2 * np.exp(- 2* a**2)
#
#zero = Np * (qt.coherent(N, a) + qt.coherent(N, -a))
#one = Np * (qt.coherent(N, 1j * a) + qt.coherent(N, -1j * a)) 
#
#a = qt.destroy(N)
#
#c1 = .2
#c2 = np.sqrt(1 - c1**2)
#
#state = c1 * zero + c2 * one
#state =  a * state
#state =  a * state
#
#print(ops.parity_measurement(N, 0, 1, state, 1))




def hypercube_state(N, vertices, a):
    theta = 2*np.pi/vertices
    state = qt.basis(N) * 0
    for i in range(vertices):
        state += qt.coherent(N, np.exp(i*1j*theta) *a)
    return state/state.norm()

N = 30
a = 8
state2 = hypercube_state(N, 3, a)
wp.plot(state2, [-10, 10])