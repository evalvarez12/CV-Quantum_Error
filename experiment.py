# -*- coding: utf-8 -*-
"""
Dummy file to test stuff

@author: Eduardo Villasenor
"""

import cv_system as cv
import qutip as qt
import numpy as np
import operations
import matplotlib.pyplot as plt
import hamiltonians as ham
import symplectic as sym
import theory

N = 2
a = qt.tensor(qt.destroy(N), qt.qeye(N))
b = qt.tensor(qt.qeye(N), qt.destroy(N))
z = np.pi/3
H_mat = ham.H_beam_splitter2(z)
print("H_mat")
print(H_mat)

H = ham.Hamiltonian(N, H_mat)

H_ref = (-1j * z * a * b.dag() + 1j * z.conjugate() * a.dag() * b)
print("H")
print(H)


S = sym.symplectic(H_mat)
S_ref = theory.beam_splitter(z)
print("-----------------")
print(S)
print(S_ref)
