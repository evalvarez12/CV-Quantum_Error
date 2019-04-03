# -*- coding: utf-8 -*-
"""
Fucntion for the different measurements of entanglement/quantumness
Created on Tue Apr  2 15:27:20 2019

@author: Eduardo Villasenor
"""

import numpy as np
import scipy.linalg as la

def key_rate(sys, f, p):
    Va = sys.get_simple_CM_V(0)
    Vb = sys.get_simple_CM_V(1)
    Ve = sys.get_simple_CM_V(2)
    Vf = sys.get_simple_CM_V(3)
    
    Cab = sys.get_simple_CM_C([0, 1])
    Cbe = sys.get_simple_CM_C([1, 2])
    Cbf = sys.get_simple_CM_C([1, 3])
    Cef = sys.get_simple_CM_C([2, 3])

    I_shared = I(Va, Vb, Cab)
    I_stolen = X(Vb, Ve, Vf, Cbe, Cbf, Cef)
    
    k_rate = p * (f * I_shared - I_stolen)
    return k_rate
    

def I(Va, Vb, Cab):
    # Conditional variance 
    Vb_a = Vb - (Cab**2)/Va
    
    # Shared information
    si = np.log2(Vb/Vb_a)/2
    return si
    
def X(Vb, Ve, Vf, Cbe, Cbf, Cef):
    I = np.eye(2)
    S = np.diag([1, -1])
    
    Mef = np.block([[I*Ve, S*Cef],[S*Cef, I*Vf]])
    vec = np.block([Cbe*I, Cbf*S])
    Mef_b = Mef - np.dot(vec.transpose(), np.dot(np.diag([1/Vb, 0]), vec))
    
    vef = symplectic_eigenvalues(Mef)
    vef_b = symplectic_eigenvalues(Mef_b)
    
    # Calculate the stolen information
    si = sum(g(vef)) - sum(g(vef_b))
    return si
    
        
def g(v):
    return ((v+1)/2)*(np.log2((v+1)/2)) - ((v-1)/2)*(np.log2((v-1)/2))    
    
    
def symplectic_eigenvalues(cm):
    N = cm.shape[0]/2
    omega = symplectic_form(N)
    eigvals = np.linalg.eigvals(1j*np.dot(omega, cm))
    return eigvals


def symplectic_form(N):
    w = np.array([[0, 1], [-1, 0]])

    # Use direct sum to calculate the symplectic form
    return la.block_diag(*[w]*N)
