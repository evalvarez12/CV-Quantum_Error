# -*- coding: utf-8 -*-
"""
Testing the different operations with symplectic

Created on Tue Apr 23 16:48:41 2019

@author: Eduardo Villasenor
"""

from hamiltonians import *
import numpy as np
import qutip as qt

N = 3

r = 2
# phase = np.pi/2*0
# z = r * np.exp(1j*phase)
z = r
theta = np.pi/6

a = qt.tensor(qt.destroy(N), qt.qeye(N))
b = qt.tensor(qt.qeye(N), qt.destroy(N))

#############################
print('-------------- TWO MODE SQUEEZING ---------------------')


H_mat = H_two_mode_squeeze(z)
H = Hamiltonian(N, H_mat)

H_ref = 1j*z.conjugate()*a*b - 1j*z*a.dag()*b.dag()
U_ref = (-1j*H_ref).expm()

U_qtref = qt.squeezing(a, b, 2*z)
print("H == H_ref:", H == H_ref)
print("U == U_ref:", U_ref == U_qtref)

S = Symplectic(H_mat)
print(np.round(S, 4))

print('Ref ------')
ch = np.cosh(-z/2)
sh = np.sinh(-z/2)
S_ref =np.array([[ch, 0, sh, 0], [0, ch, 0, -sh], [sh, 0, ch, 0], [0, -sh, 0, ch]])
print(S_ref)

print('--------------- SINGLE MODE SQUEEZING --------------------')

H_mat = H_single_mode_squeeze(z)
H = Hamiltonian(N, H_mat)

H_ref = (1j*z.conjugate()*a**2 - 1j*z*a.dag()**2)/2
U_ref = (-1j*H_ref).expm()

U_qtref = qt.squeezing(a, a, z)
print("H == H_ref:", H == H_ref)
print("U == U_ref:", U_ref == U_qtref)

S = Symplectic(H_mat)
print(np.round(S, 4))

print('Ref ------')
phasee = np.pi/2*0
cc = np.cos(phasee)
ss = np.sin(phasee)
ch = np.cosh(-r/2)
sh = np.sinh(-r/2)
S_ref =np.array([[ch + cc*sh, ss*sh, 0, 0], [ss*sh, ch- cc*sh, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
print(S_ref)

print('--------------- PHASE SHIFT --------------------')

H_mat = H_phase_shift(theta)
H = Hamiltonian(N, H_mat)

H_ref = -theta*(a * a.dag() + a.dag() * a)
U_ref = (-1j*H_ref).expm()

# U_qtref = qt.phase_gate(N, theta)
print("H == H_ref:", H == H_ref)
# print(U_ref == U_qtref)

S = Symplectic(H_mat)
print(np.round(S, 4))

print('Ref ------')
cc = np.cos(theta)
ss = np.sin(theta)
S_ref =np.array([[cc, -ss, 0, 0], [ss, cc, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
print(S_ref)

Sps = S_ref
print('--------------- BEAM SPLITTER --------------------')
z = np.pi/4

H_mat = H_beam_splitter(z)
H = Hamiltonian(N, H_mat)

H_ref = (1j*z.conjugate()*a * b.dag() -1j*z* a.dag() * b + 1j*z.conjugate()*b.dag()*a - 1j*z*b*a.dag())
# H_ref = -theta*(a * b.dag() + b.dag() * a)

U_ref = (-1j*H_ref).expm()

# U_qtref = qt.phase_gate(N, theta)
print("H == H_ref:", H == H_ref)
# print(U_ref == U_qtref)

S = Symplectic(H_mat)
print(np.round(S, 4))

print('Ref ------')
c1 = np.cos(np.absolute(z))
s1 = np.sin(np.absolute(z))
c2 = np.cos(np.angle(z))
s2 = np.sin(np.angle(z))

S_ref =np.array([[c1, 0, s1*c2, s1*s2], [0, c1, -s1*s2, s1*c2], [-s1*c2, -s1*s2, c1, 0], [s1*s2, -s1*c2, 0, c1]])
print(S_ref)
