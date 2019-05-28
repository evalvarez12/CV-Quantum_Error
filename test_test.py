# -*- coding: utf-8 -*-
"""
Dummy file to test stuff

@author: Eduardo Villasenor
"""

import cv_system as cv
import qutip as qt
import numpy as np
import operations
# Define some initial parameters
N = 30
t = .3
theta = np.arccos(np.sqrt(t))

alpha = 2 + 2j
beta = 3 + 1j

a = qt.tensor(qt.destroy(N), qt.identity(N))
b = qt.tensor(qt.identity(N), qt.destroy(N))

state_a = (alpha*a.dag() - alpha.conjugate()*a).expm()
state_b = (beta*b.dag() - beta.conjugate()*b).expm()

U = operations.beam_splitter(N, theta)


a_out = np.cos(theta)*a + -1*np.sin(theta)*b
b_out = np.cos(theta)*b + -1*np.sin(theta)*a

state_a_ref = (alpha*a_out.dag() - alpha.conjugate()*a_out).expm()
state_b_ref = (beta*b_out.dag() - beta.conjugate()*b_out).expm()

state_a_out = U * state_a * U.dag()
state_b_out = U * state_b * U.dag()

print(state_a_out == state_a_ref)
print(state_b_out == state_b_ref)

#N = 10
#sys = cv.System(N, 1)
#mu = 1.1
#r = np.arcsinh(np.sqrt(mu))
#sys.apply_SMD(r)
#
#ref_state = sys.state
#k = .1
#m_aux = .01
#r_aux = np.arcsinh(np.sqrt(m_aux))
##p = sys.apply_scissors(k, r_aux)
#p = sys.apply_scissors_options(k, r_aux, 0, 'c')
#
##sys.apply_photon_subtraction(k, 0)
#g = np.sqrt((1-k)/k)
#p_ref = np.exp((1-g**2)*np.linalg.norm(r)**2)/(1 + g**2)
#
#
#print(ref_state)
#print(sys.state)
#
#print("g:", g)
#print("p:", p)
#print("p_ref:", p_ref)
