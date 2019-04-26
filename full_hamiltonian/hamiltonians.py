# -*- coding: utf-8 -*-
"""
Symplectic transformations for common operations

Created on Tue Apr 23 10:56:50 2019

@author: Eduardo Villasenor
"""

import numpy as np
import qutip as qt


def beam_spliter(N, z):
    return Hamiltonian(N, H_beam_splitter(z))


def phase_shift(N, theta):
    return Hamiltonian(N, H_phase_shift(theta))


def sigle_mode_squeeze(N, z):
    return Hamiltonian(N, H_single_mode_squeeze(z))


def two_mode_squeeze(N, z):
    return Hamiltonian(N, H_two_mode_squeeze(z))


def H_single_mode_squeeze(z):
    H = np.zeros([4, 4], dtype=complex)
    H[2, 0] = (-1j * z.conjugate())/2
    H[0, 2] = (1j * z)/2
    return H


def H_two_mode_squeeze(z):
    H = np.zeros([4, 4], dtype=complex)
    H[3, 0] = (-1j * z.conjugate())/2
    H[0, 3] = (1j * z)/2

    H[2, 1] = (-1j * z.conjugate())/2
    H[1, 2] = (1j * z)/2
    return H


def H_beam_splitter(z):
    H = np.zeros([4, 4], dtype=complex)
    H[1, 0] = (-1j * z.conjugate())/2
    H[0, 1] = (1j * z)/2

    H[3, 2] = (1j * z)/2
    H[2, 3] = (-1j*z.conjugate())/2
    return H


def H_phase_shift(theta):
    H = np.zeros([4, 4], dtype=complex)
    H[0, 0] = -theta/2
    H[2, 2] = -theta/2
    return H


def Hamiltonian(N, H_mat):
    a = qt.tensor(qt.destroy(N), qt.qeye(N))
    b = qt.tensor(qt.qeye(N), qt.destroy(N))

    vec = np.array([a, b, a.dag(), b.dag()]).transpose()
    vec_dag = np.array([a.dag(), b.dag(), a, b])
    H = np.dot(vec_dag, np.dot(H_mat, vec))
    return H

