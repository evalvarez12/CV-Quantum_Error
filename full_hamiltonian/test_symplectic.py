# -*- coding: utf-8 -*-
"""
Testing the different operations with symplectic

Created on Tue Apr 23 16:48:41 2019

@author: Eduardo Villasenor
"""

from symplectic import *
import numpy as np
import qutip as qt

N = 3

r = 2
phase = np.pi/2*0
z = r * np.exp(1j*phase)
z = r

a = qt.tensor(qt.destroy(N), qt.qeye(N))
b = qt.tensor(qt.qeye(N), qt.destroy(N))

#############################


H_mat = H_two_mode_squeeze(z * np.exp(-1j*np.pi))
H = Hamiltonian(N, H_mat)

H_ref = 1j*z.conjugate()*a*b - 1j*z*a.dag()*b.dag()
U_ref = (-1j*H_ref).expm()

U_qtref = qt.squeezing(a, b, 2*z)
print(H == H_ref)
print(U_ref == U_qtref)

S = Symplectic(H_mat)
print(S)

print('Ref ------')
ch = np.cosh(z)
sh = np.sinh(z)
S_ref =np.array([[ch, 0, sh, 0], [0, ch, 0, -sh], [sh, 0, ch, 0], [0, -sh, 0, ch]])
print(S_ref)

print('-----------------------------------')

H_mat = H_single_mode_squeeze(z)
H = Hamiltonian(N, H_mat)

H_ref = (1j*z.conjugate()*a**2 - 1j*z*a.dag()**2)
U_ref = (-1j*H_ref).expm()

U_qtref = qt.squeezing(a, a, 2*z)
print(H == H_ref)
print(U_ref == U_qtref)

S = Symplectic(H_mat)
print(S)

print('Ref ------')
phasee = np.pi/2*0
cc = np.cos(phasee)
ss = np.sin(phasee)
ch = np.cosh(-r)
sh = np.sinh(-r)
S_ref =np.array([[ch + cc*sh, ss*sh, 0, 0], [ss*sh, ch- cc*sh, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
print(S_ref)


