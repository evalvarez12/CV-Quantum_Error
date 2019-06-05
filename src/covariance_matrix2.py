# -*- coding: utf-8 -*-
"""
Functions to work on CV systems using covariance matrices

Note  r = (x1, x2, ...,xN, p1, p2, ..., pN) and direct sumation is used
though the code


Created on Wed Mar 13 11:04:05 2019

@author: Eduardo Villasenor
"""


import numpy as np
import tools_cm as ts

#
#a = np.array([[1,2],[2,3]])
#b = a*0
#np.block([[a,b],[b,a]])


def beam_splitter_symplectic(t, positions=[0,1], nmodes=2):
    # The symplectic transformation for  a single mode
    S = np.array([[np.sqrt(t), np.sqrt(1 - t)], [-np.sqrt(1 - t), np.sqrt(t)]])
    
    # Apply the direct sum to obtain the full transformation
    S = ts.direct_sum_singles([S, S], positions, nmodes)
    return S


def pure_loss(cm, pos, nblocks, eta):
    X = np.sqrt(eta) * np.identity(2)
    Y = (1 - eta) * np.identity(2)
    
    X = ts.direct_sum_singles([X], [pos], nblocks)
    Y = ts.direct_sum_singles([Y], [pos], nblocks)
    
    return np.dot(np.transpose(X), np.dot(cm, X)) + Y



def tmsv_cm(s):
    # TMSV states are always subsequent in their indices
    Vp = np.array([[np.cosh(2*s), np.sinh(2*s)], [np.sinh(2*s), np.cosh(2*s)]])
    Vm = np.array([[np.cosh(2*s), -np.sinh(2*s)], [-np.sinh(2*s), np.cosh(2*s)]])
    zeros = np.zeros((2, 2))
    
    return np.block([[Vp, zeros], [zeros, Vm]])


def measurement(cm, mode, Nmodes, cm_m, r_m):
    a, b = modes
    
    # Grab the rest of the CM
    ind_start = (mode - 1) * N
    ind_end = ind_start + N
    cm_block = cm[ind_start:ind_end, ind_start:ind_end]
    
    
        
def measurement_probabilty(cm, s, cm_m, r_m):
    a = np.linalg.inv(cm + cm_m)
    a = np.dot(r_m-s , np.dot(a,(r_m - s)))
    p = np.exp(-a)/np.sqrt(np.linalg.det(cm + cm_m))
    
    return p 


def symplectic_eigenvalues(cm, Nmodes):
    omega = symplectic_form(Nmodes)
    eigvals = np.linalg.eigvals(1j*np.dot(omega, cm))
    return eigvals



def symplectic_form(Nmodes):
    w = np.array([[0, 1], [-1, 0]])
    
    # Use direct sum to calculate the symplectic form
    return ts.direct_sum([w]*Nmodes)



s = 10 + 1j

state = tmsv_cm(s)
eigvals = symplectic_eigenvalues(state, 2)

print(eigvals)





#a = np.array([[1,2],[3,4]])
#
#print(direct_sum([a,a], [0,3], 3))
#
#print(beam_splitter_symplectic(np.pi/4, [0, 1], 2))
    
#w = np.array([[0, 1], [-1, 0]])
#
#print(direct_sum([w, w, w], [0, 1, 2], 3))