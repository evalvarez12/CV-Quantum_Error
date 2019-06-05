# -*- coding: utf-8 -*-
"""
Symplectic representation of evolution operators
Created on Wed Apr 24 14:54:56 2019

@author: Eduardo Villasenor 
"""

import numpy as np
import scipy.linalg as la
from . import hamiltonians as ham
from . import tools



def beam_splitter(z, pos=[0, 1], Nmodes=2):
    H = ham.H_beam_splitter(z)
    S = symplectic(H)
    
    if Nmodes > 2:
        S = tools.reorder_two_mode_symplectic(S, pos, Nmodes)
    return S
    

def phase_shift(z, pos=0, Nmodes=1):
    H = ham.H_phase_shift(z)
    S = symplectic(H)
    
    # Keep only the first mode on which the phase shift acts
    S = S[0:2, 0:2]
    
    if Nmodes > 1:
        S = tools.direct_sum_singles([S], [pos], Nmodes)
    return S


def single_mode_squeeze(z, pos=0, Nmodes=1):
    H = ham.H_single_mode_squeeze(z)
    # Extra factor of 2 requied specifically for this case
    S = symplectic(H)

    # Keep only the first mode on which the squeezing acts
    S = S[0:2, 0:2]
    
    if Nmodes > 1:
        S = tools.direct_sum_singles([S], [pos], Nmodes)
    return S


def two_mode_squeeze(z, pos=[0,1], Nmodes=2):
    H = ham.H_two_mode_squeeze(z)
    S = symplectic(H)
    
    if Nmodes > 2:
        S = tools.reorder_two_mode_symplectic(S, pos, Nmodes)
    return S


def symplectic(H_mat):
    K = np.block([[np.eye(2), np.zeros([2, 2])], [np.zeros([2, 2]), -np.eye(2)]])
    # Extra factor on the Hamiltonian matrix requied to get correct representation
    S = la.expm(-1j*np.dot(K, 2 * H_mat))

#    size = H_mat.shape[0]
    L = quad_basis_transform()
    T = quad_basis_reorder()
    S = np.dot(L.transpose().conjugate(), np.dot(S, L))
#    return S
    return np.dot(T, np.dot(S, T)).real


def quad_basis_transform():
#    if modes == 4:
    L = np.block([[np.eye(2), 1j*np.eye(2)], [np.eye(2), -1j*np.eye(2)]])/np.sqrt(2)
#    else:
#    L = np.array([[1, 1j], [1, -1j]])/np.sqrt(2)
    return L


def quad_basis_reorder():
    T = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    return T
