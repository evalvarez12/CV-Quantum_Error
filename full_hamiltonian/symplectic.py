# -*- coding: utf-8 -*-
"""
Symplectic transformations for common operations

Created on Tue Apr 23 10:56:50 2019

@author: Eduardo Villasenor
"""

import numpy as np
import qutip as qt
import scipy.linalg as la


def H_single_mode_squeeze(z):
    H = np.zeros([4, 4], dtype=complex)
    H[2, 0] = (1j * z.conjugate())/2
    H[0, 2] = (-1j * z)/2

#    H[2, 0] = z.conjugate()
#    H[0, 2] = z
    return H


def H_two_mode_squeeze(z):
    H = np.zeros([4, 4], dtype=complex)
    H[3, 0] = (1j * z.conjugate())/2
    H[0, 3] = (-1j * z)/2

    H[2, 1] = (1j * z.conjugate())/2
    H[1, 2] = (-1j * z)/2
    return H


def H_beam_splitter(theta):
    H = np.zeros([4, 4])
    H[1, 0] = -theta
    H[0, 1] = -theta

    H[3, 2] = -theta
    H[2, 3] = -theta
    return H


def H_phase_shift(theta):
    H = np.zeros([4, 4])
    H[0, 0] = -theta
    H[2, 2] = -theta
    return H


def Hamiltonian(N, H_mat):
    a = qt.tensor(qt.destroy(N), qt.qeye(N))
    b = qt.tensor(qt.qeye(N), qt.destroy(N))

    vec = np.array([a, b, a.dag(), b.dag()]).transpose()
    vec_dag = np.array([a.dag(), b.dag(), a, b])
    H = np.dot(vec_dag, np.dot(H_mat, vec))
    return H


def Symplectic(H_mat):
    K = np.block([[np.eye(2), np.zeros([2, 2])], [np.zeros([2, 2]), -np.eye(2)]])
    S = la.expm(-1j*np.dot(K, H_mat))

#    size = H_mat.shape[0]
    L = quad_basis_transform()
    T = quad_basis_reorder()
    S = np.dot(L.transpose().conjugate(), np.dot(S, L))
#    return S
    return np.dot(T, np.dot(S, T))


def quad_basis_transform():
#    if modes == 4:
    L = np.block([[np.eye(2), 1j*np.eye(2)], [np.eye(2), -1j*np.eye(2)]])/np.sqrt(2)
#    else:
#    L = np.array([[1, 1j], [1, -1j]])/np.sqrt(2)
    return L


def quad_basis_reorder():
    T = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    return T
