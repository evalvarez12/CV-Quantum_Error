# -*- coding: utf-8 -*-
"""
Testing the different operations with symplectic

Created on Tue Apr 23 16:48:41 2019

@author: Eduardo Villasenor
"""

import hamiltonians as ham
import symplectic as sym
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
        H_mat = ham.H_two_mode_squeeze(z)
        H = ham.Hamiltonian(N, H_mat)


        H_ref = -1j * z.conjugate() * a * b  + 1j * z * a.dag() * b.dag()
        U_ref = (-1j * H_ref).expm()

        # Qutip needs a -1 and in TMS a factor of 2 in squeezing parameter
        U_qtref = qt.squeezing(a, b, -2*z)

        S = sym.symplectic(H_mat)
        S_ref = theory.two_mode_squeeze(z)

        self.assertTrue(U_ref == U_qtref)
        self.assertTrue(H == H_ref)
        np.testing.assert_array_almost_equal(S, S_ref)


    def test_single_mode_squeezing(self):
        z = 2 * np.exp(1j*np.pi/6)
        H_mat = ham.H_single_mode_squeeze(z)
        H = ham.Hamiltonian(N, H_mat)

        H_ref = (-1j * z.conjugate() * a**2 + 1j * z * a.dag()**2)/2
        U_ref = (-1j * H_ref).expm()

        # Qutip needs a -1 in squeezing paramenter
        U_qtref = qt.squeezing(a, a, -z)

        S = sym.symplectic(H_mat)
        S_ref = theory.single_mode_squeeze_2modes(z)

        self.assertTrue(U_ref == U_qtref)
        self.assertTrue(H == H_ref)
        np.testing.assert_array_almost_equal(S, S_ref)



    def test_phase_shift(self):
        theta = np.pi/6
        H_mat = ham.H_phase_shift(theta)
        H = ham.Hamiltonian(N, H_mat)

        H_ref = -theta*(a * a.dag() + a.dag() * a)/2

        S = sym.symplectic(H_mat)
        S_ref = theory.phase_shift_2modes(theta)

        self.assertTrue(H == H_ref)
        np.testing.assert_array_almost_equal(S, S_ref)


    def test_beam_splitter(self):
        z = np.pi/4 * np.exp(1j*np.pi/4)
        H_mat = ham.H_beam_splitter(z)
        H = ham.Hamiltonian(N, H_mat)

        H_ref = (-1j * z * a * b.dag() + 1j * z.conjugate() * a.dag() * b)

        S = sym.symplectic(H_mat)
        S_ref = theory.beam_splitter(z)

        self.assertTrue(H == H_ref)
        np.testing.assert_array_almost_equal(S, S_ref, decimal=5)


    def test_single_mode_squeezing_S(self):
         z = 2 * np.exp(1j*np.pi/6)
         S = sym.single_mode_squeeze(z)
         S_ref = theory.single_mode_squeeze(z)
         np.testing.assert_array_almost_equal(S, S_ref, decimal=5)
         

    def test_two_mode_squeezing_S(self):
         z = 2 * np.exp(1j*np.pi/6)
         S = sym.two_mode_squeeze(z)
         S_ref = theory.two_mode_squeeze(z)
         np.testing.assert_array_almost_equal(S, S_ref, decimal=5)
         

    def test_phase_shift_S(self):
         theta = np.pi/6
         S = sym.phase_shift(theta)
         S_ref = theory.phase_shift(theta)
         np.testing.assert_array_almost_equal(S, S_ref, decimal=5)
         
         
    def test_beam_splitter_S(self):
         z = 2 * np.exp(1j*np.pi/6)
         S = sym.beam_splitter(z)
         S_ref = theory.beam_splitter(z)
         np.testing.assert_array_almost_equal(S, S_ref, decimal=5)

if __name__ == '__main__':
    unittest.main()
