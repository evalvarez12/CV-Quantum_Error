# -*- coding: utf-8 -*-
"""
Fucntion for the different measurements of entanglement/quantumness
Created on Tue Apr  2 15:27:20 2019

@author: Eduardo Villasenor
"""

import numpy as np
import scipy.linalg as la

def key_rate(sys, f, p):
    Va = sys.get_simple_CM_V(0).norm()
    Vb = sys.get_simple_CM_V(1).norm()
    Ve = sys.get_simple_CM_V(2).norm()
    Vf = sys.get_simple_CM_V(3).norm()
    
    Cab = sys.get_simple_CM_C([0, 1]).norm()
    Cbe = sys.get_simple_CM_C([1, 2]).norm()
    Cbf = sys.get_simple_CM_C([1, 3]).norm()
    Cef = sys.get_simple_CM_C([2, 3]).norm()

    I_shared = I(Va, Vb, Cab)
    I_stolen = X(Vb, Ve, Vf, Cbe, Cbf, Cef)
    
    k_rate = p * (f * I_shared - I_stolen)
    return k_rate
    

def key_rate_compare(sys, f, p, mpnA, mpnE, t):
    Va = sys.get_simple_CM_V(0)
    Vb = sys.get_simple_CM_V(1)
    Ve = sys.get_simple_CM_V(2)
    Vf = sys.get_simple_CM_V(3)
    
    Cab = sys.get_simple_CM_C([0, 1])
    Cbe = sys.get_simple_CM_C([1, 2])
    Cbf = sys.get_simple_CM_C([1, 3])
    Cef = sys.get_simple_CM_C([2, 3])

    # NOPS
    Va2 = 2 * mpnA + 1
    Vf2 = 2 * mpnE + 1
    Vb2 = t * Va2 + (1 - t) * Vf2
    Ve2 = (1 - t) * Va2 + t * Vf2
    
    Cab2 = np.sqrt(t) * 2 * np.sqrt(mpnA**2 + mpnA)
    Cef2 = np.sqrt(t) * 2 * np.sqrt(mpnE**2 + mpnE)
    Cbf2 = np.sqrt(1 - t)* 2 * np.sqrt(mpnE**2 + mpnE) 
    Cbe2 = np.sqrt(t * (1 - t)) * (Vf2 - Va2) 
    
    # TPS
#    Va = 1 + 2 * (N * mpnA + N + mpnA * Tbs)/(1 + mpnA - mpnA * Tbs);
#    Vb1 = 1 + 2 * (N + 1) * mpnA * Tbs / (1 + mpnA - mpnA * Tbs);
#    Vf = 1 + 2 * mpnE;
#    Vf = t * Vb1 + (1 - t) * Vf;
#    Ve = (1 - t) * Vb1 + t * Vf;
#
#    Cab = np.sqrt(t) * 2 * np.sqrt(mpnA * Tbs/(1 + mpnA)) * (N + 1) * (1 + mpnA)/(1 + mpnA - mpnA * Tbs);
#    Cef = np.sqrt(t) * 2 * np.sqrt(mpnE^2 + mpnE);
#    Cfb = np.sqrt(1 - t) * 2 * np.sqrt(mpnE^2 + mpnE);
#    Ceb = np.sqrt(t* (1 - t)) * (Vf - Vb1);

    sys.set_quadratures_basis() 
    Va3 = sys.get_CM_entry([0, 0])
    Vb3 = sys.get_CM_entry([2, 2])
    Ve3 = sys.get_CM_entry([4, 4])
    Vf3 = sys.get_CM_entry([6, 6])
    
    Cab3 = sys.get_CM_entry([0, 2])
    Cef3 = sys.get_CM_entry([4, 6])
    Cbf3 = sys.get_CM_entry([2, 7])
    Cbe3 = sys.get_CM_entry([2, 5])


#    Cab2 = abs(Cab2)
#    Cef2 = abs(Cef2)
#    Cbf2 = abs(Cbf2)
#    Cbe2 = abs(Cbe2)

    print("----------------t---------------", t)
    print("Va:", Va, Va2, Va3)
    print("Vb:", Vb, Vb2, Vb3)
    print("Ve:", Ve, Ve2, Ve3)
    print("Vf:", Vf, Vf2, Vf3)
    print("Cab:", Cab, Cab2, Cab3)
    print("Cbe:", Cbe, Cbe2, Cbe3)
    print("Cbf:", Cbf, Cbf2, Cbf3)
    print("Cef:", Cef, Cef2, Cef3)

#    CM = sys.get_full_CM()
#    print("CM:", CM)

    I_shared = I(Va3, Vb3, Cab3)
    I_stolen = X(Vb3, Ve3, Vf3, Cbe3, Cbf3, Cef3)
    
    k_rate = p * (f * I_shared - I_stolen)
    
#    print("I_shared:", f*I_shared)
#    print("I_stolen:", I_stolen)
    
    return k_rate

def key_rate_NOsimple(mpnA, mpnE, t):
    Va = 2 * mpnA + 1
    Vf = 2 * mpnE + 1
    Vb = t * Va + (1 - t) * Vf
    Ve = (1 - t) * Va + t * Vf
    
    Cab = np.sqrt(t) * 2 * np.sqrt(mpnA**2 + mpnA)
    Cef = np.sqrt(t) * 2 * np.sqrt(mpnE**2 + mpnE)
    Cbf = np.sqrt(1 -t)* 2 * np.sqrt(mpnE**2 + mpnE)
    Cbe = np.sqrt(t * (1 - t)) * (Vf - Va)
    
    I_shared = I(Va, Vb, Cab)
    I_stolen = X(Vb, Ve, Vf, Cbe, Cbf, Cef)
    
    p = 1
    f = 1
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
#    vef2 = .5 * np.sqrt((Ve + Vf)**2 - 4 * Cef**2) + .5 * np.array([Ve-Vf, Vf -Ve])
#    print("Vef:", vef, vef2)
    
    vef_b = symplectic_eigenvalues(Mef_b)
    
    # Calculate the stolen information
    si = sum(g(vef)) - sum(g(vef_b))
    si = si.real/2
    return si
    
        
def g(v):
    return ((v+1)/2)*(np.log2((v+1)/2)) - ((v-1)/2)*(np.log2((v-1)/2))    
    
    
def symplectic_eigenvalues(cm):
    N = int(cm.shape[0]/2)
    omega = symplectic_form(N)
    eigvals = np.linalg.eigvals(1j*np.dot(omega, cm))
    return eigvals


def symplectic_form(N):
    w = np.array([[0, 1], [-1, 0]])

    # Use direct sum to calculate the symplectic form
    return la.block_diag(*[w]*N)
