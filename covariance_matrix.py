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


def beam_splitter_symplectic(theta, positions=[0,1], nmodes=2):
    # The symplectic transformation for  a single mode
    S = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    
    # Apply the direct sum to obtain the full transformation
    S = direct_sum([S, S], positions, nmodes)
    return S


def tmsv(s):
    # TMSV states are always subsequent in their indices
    Vp = np.array([[np.cosh(2*s), np.sinh(2*s)], [np.sinh(2*s), np.cosh(2*s)]])
    Vm = np.array([[np.cosh(2*s), -np.sinh(2*s)], [-np.sinh(2*s), np.cosh(2*s)]])
    zeros = np.zeros((2, 2))
    
    return np.block([[Vp, zeros], [zeros, Vm]])


def measurement(cm, modes, N, cm_m, r_m):
    a, b = modes
    cm_block = cm[a:b, a:b]
    
    
    
    






#a = np.array([[1,2],[3,4]])
#
#print(direct_sum([a,a], [0,3], 3))
#
#print(beam_splitter_symplectic(np.pi/4, [0, 1], 2))
    
#w = np.array([[0, 1], [-1, 0]])
#
#print(direct_sum([w, w, w], [0, 1, 2], 3))