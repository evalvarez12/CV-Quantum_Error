# -*- coding: utf-8 -*-
"""
Testing the different operations with symplectic

Created on Tue Apr 23 16:48:41 2019

@author: Eduardo Villasenor
"""

from hamiltonians import *
from symplectic import *
import theory
import numpy as np
import qutip as qt
import unittest


N = 10
a = qt.tensor(qt.destroy(N), qt.qeye(N))
b = qt.tensor(qt.qeye(N), qt.destroy(N))

class TestSymplecticMethods(unittest.TestCase):

    def test_two_mode_squeezing(self):
        z = 2 * np.exp(1j*np.pi/6)
        H_mat = H_two_mode_squeeze(z)
        H = Hamiltonian(N, H_mat)


        H_ref = -1j * z.conjugate() * a * b  + 1j * z * a.dag() * b.dag()
        U_ref = (-1j * H_ref).expm()

        # Qutip needs a -1 and in TMS a factor of 2 in squeezing parameter
        U_qtref = qt.squeezing(a, b, -2*z)

        S = symplectic(H_mat)
        S_ref = theory.two_mode_squeeze(z)

#        print(np.round(S, 5))
#        print(np.round(S_ref, 5))
        self.assertTrue(U_ref == U_qtref)
        self.assertTrue(H == H_ref)
        np.testing.assert_array_almost_equal(S, S_ref)


    def test_single_mode_squeezing(self):
        z = 2 * np.exp(1j*np.pi/6)
        H_mat = H_single_mode_squeeze(z)
        H = Hamiltonian(N, H_mat)

        H_ref = (-1j * z.conjugate() * a**2 + 1j * z * a.dag()**2)/2
        U_ref = (-1j * H_ref).expm()

        # Qutip need a -1 in squeezing paramenter
        U_qtref = qt.squeezing(a, a, -z)

        S = symplectic(H_mat)
        S_ref = theory.single_mode_squeeze(z)

#        print(np.round(S, 5))
#        print(np.round(S_ref, 5))
        self.assertTrue(U_ref == U_qtref)
        self.assertTrue(H == H_ref)
        np.testing.assert_array_almost_equal(S, S_ref)



    def test_phase_shift(self):
        theta = np.pi/6
        H_mat = H_phase_shift(theta)
        H = Hamiltonian(N, H_mat)

        H_ref = -theta*(a * a.dag() + a.dag() * a)/2

        S = symplectic(H_mat)
        S_ref = theory.phase_shift(theta)

#        print(np.round(S, 5))
#        print(np.round(S_ref, 5))
        self.assertTrue(H == H_ref)
        np.testing.assert_array_almost_equal(S, S_ref)


    def test_beam_splitter(self):
        z = np.pi/4 * np.exp(1j*np.pi/4)
        H_mat = H_beam_splitter(z)
        H = Hamiltonian(N, H_mat)

        H_ref = (-1j * z * a * b.dag() + 1j * z.conjugate() * a.dag() * b)

        S = symplectic(H_mat)
        S_ref = theory.beam_splitter(z)

#        print(np.round(S, 5))
#        print(np.round(S_ref, 5))
        self.assertTrue(H == H_ref)
        np.testing.assert_array_almost_equal(S, S_ref, decimal=5)


if __name__ == '__main__':
    unittest.main()

#############################
# print('-------------- TWO MODE SQUEEZING ---------------------')
#
# z = 2 * np.exp(1j*np.pi/4)
# H_mat = H_two_mode_squeeze(z)
# H = Hamiltonian(N, H_mat)
#
#
# H_ref = -1j * z.conjugate() * a * b  + 1j * z * a.dag() * b.dag()
# U_ref = (-1j * H_ref).expm()
#
# # Qutip needs a -1 and in TMS a factor of 2 in squeezing parameter
# U_qtref = qt.squeezing(a, b, -2*z)
# print("H == H_ref:", H == H_ref)
# print("U2 == U_ref:", U_ref == U_qtref)
#
# S = symplectic(H_mat)
# print(np.round(S, 5))
#
# print('Ref ------')
# S_ref = theory.two_mode_squeeze(z)
# print(np.round(S_ref, 5))
#
#
# print('--------------- SINGLE MODE SQUEEZING --------------------')
#
# H_mat = H_single_mode_squeeze(z)
# H = Hamiltonian(N, H_mat)
#
# H_ref = (-1j * z.conjugate() * a**2 + 1j * z * a.dag()**2)/2
# U_ref = (-1j * H_ref).expm()
#
# # Qutip need a -1 in squeezing paramenter
# U_qtref = qt.squeezing(a, a, -z)
# print("H == H_ref:", H == H_ref)
# print("U == U_ref:", U_ref == U_qtref)
#
# S = symplectic(H_mat)
# print(np.round(S, 4))
#
# print('Ref ------')
# S_ref = theory.single_mode_squeeze(z)
# print(np.round(S_ref, 5))
#
#
# print('--------------- PHASE SHIFT --------------------')
#
# theta = np.pi/6
# H_mat = H_phase_shift(theta)
# H = Hamiltonian(N, H_mat)
#
# H_ref = -theta*(a * a.dag() + a.dag() * a)/2
# U_ref = (-1j * H_ref).expm()
#
# # U_qtref = qt.phase_gate(N, theta)
# print("H == H_ref:", H == H_ref)
# # print(U_ref == U_qtref)
#
# S = symplectic(H_mat)
# print(np.round(S, 4))
#
# print('Ref ------')
# S_ref = theory.phase_shift(theta)
# print(np.round(S_ref, 5))
#
# Sps = S_ref
#
#
# print('--------------- BEAM SPLITTER --------------------')
# z = np.pi/4
#
# H_mat = H_beam_splitter(z)
# H = Hamiltonian(N, H_mat)
#
# H_ref = (-1j * z.conjugate() * a * b.dag() + 1j * z * a.dag() * b)
# U_ref = (-1j * H_ref).expm()
#
# # U_qtref = qt.phase_gate(N, theta)
# print("H == H_ref:", H == H_ref)
# # print(U_ref == U_qtref)
#
# S = symplectic(H_mat)
# print(np.round(S, 4))
#
# print('Ref ------')
# S_ref = theory.beam_splitter(z)
# print(np.round(S_ref, 5))
