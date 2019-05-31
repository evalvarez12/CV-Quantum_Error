# -*- coding: utf-8 -*-
"""
Fucntion for the different measurements of entanglement/quantumness
Created on Tue Apr  2 15:27:20 2019

@author: Eduardo Villasenor
"""

import numpy as np
import scipy.linalg as la
import qutip as qt

def key_rate_simple(sys, f, p):
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
    

def key_rate(sys, f, p):
    sys.set_quadratures_basis() 
    Va = sys.get_CM_entry([0, 0])
    Vb = sys.get_CM_entry([2, 2])
    Ve = sys.get_CM_entry([4, 4])
    Vf = sys.get_CM_entry([6, 6])
    
    Cab = sys.get_CM_entry([0, 2])
    Cbe = sys.get_CM_entry([2, 4])
    Cbf = sys.get_CM_entry([2, 6])
    Cef = sys.get_CM_entry([4, 6])

    I_shared = I(Va, Vb, Cab)
    I_stolen = X(Vb, Ve, Vf, Cbe, Cbf, Cef)
    
    k_rate = p * (f * I_shared - I_stolen)
    print("I_shared:", I_shared)
    print("I_stolen:", I_stolen)
    return k_rate

def key_rate_compare(sys, f, p, mpnA, mpnE, t):
    Va = sys.get_simple_CM_V(0)
    Vb = sys.get_simple_CM_V(1)
    Ve = sys.get_simple_CM_V(2)
    Vf = sys.get_simple_CM_V(3)
    
    Cab = sys.get_simple_CM_C([0, 1])
    Cbe = sys.get_simple_CM_C([2, 1])
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
#    N = sys.N
#    Va2 = 1 + 2 * (N * mpnA + N + mpnA * t)/(1 + mpnA - mpnA * t)
#    Vb2 = 1 + 2 * (N + 1) * mpnA * t / (1 + mpnA - mpnA * t)
#    Vf2 = t * Vb2 + (1 - t) * (1 + 2 * mpnE)
#    Ve2 = (1 - t) * Vb2 + t * Vf2;
#
#    Cab2 = np.sqrt(t) * 2 * np.sqrt(mpnA * t/(1 + mpnA)) * (N + 1) * (1 + mpnA)/(1 + mpnA - mpnA * t);
#    Cef2 = np.sqrt(t) * 2 * np.sqrt(mpnE**2 + mpnE);
#    Cbf2 = np.sqrt(1 - t) * 2 * np.sqrt(mpnE**2 + mpnE);
#    Cbe2 = np.sqrt(t* (1 - t)) * (Vf2 - Vb2);

    sys.set_quadratures_basis() 
    Va3 = sys.get_CM_entry([0, 0])
    Vb3 = sys.get_CM_entry([2, 2])
    Ve3 = sys.get_CM_entry([4, 4])
    Vf3 = sys.get_CM_entry([6, 6])
    
    Cab3 = sys.get_CM_entry([0, 2])
    Cbe3 = sys.get_CM_entry([2, 4])
    Cbf3 = sys.get_CM_entry([2, 6])
    Cef3 = sys.get_CM_entry([4, 6])


#    Cab2 = abs(Cab2)
#    Cef2 = abs(Cef2)
#    Cbf2 = abs(Cbf2)
#    Cbe2 = abs(Cbe2)
    
    print("----------------t---------------", t)
    print("Parameter - Simple - Theory NOPS - Calculated")
    print("Va:", Va, Va2, Va3)
    print("Vb:", Vb, Vb2, Vb3)
    print("Ve:", Ve, Ve2, Ve3)
    print("Vf:", Vf, Vf2, Vf3)
    print("Cab:", Cab, Cab2, Cab3)
    print("Cbe:", Cbe, Cbe2, Cbe3)
    print("Cbf:", Cbf, Cbf2, Cbf3)
    print("Cef:", Cef, Cef2, Cef3)

#    print("----------------t---------------", t)
#    print("Parameter - Simple - Calculated")
#    print("Va:", Va, Va3)
#    print("Vb:", Vb, Vb3)
#    print("Ve:", Ve, Ve3)
#    print("Vf:", Vf, Vf3)
#    print("Cab:", Cab, Cab3)
#    print("Cbe:", Cbe, Cbe3)
#    print("Cbf:", Cbf, Cbf3)
#    print("Cef:", Cef, Cef3)

#    CM = sys.get_full_CM()
#    print("CM:", CM)

    I_shared = I(Va3, Vb3, Cab3)
    I_stolen = X(Vb3, Ve3, Vf3, Cbe3, Cbf3, Cef3)

#    I_shared = I(Va, Vb, Cab)
#    I_stolen = X(Vb, Ve, Vf, Cbe, Cbf, Cef)
    
    k_rate = p * (f * I_shared - I_stolen)
    
#    print("I_shared:", f*I_shared)
#    print("I_stolen:", I_stolen)
    
    return k_rate

def key_rate_theory(mpnA, mpnE, t):
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
    si = np.sum(g(vef)) - np.sum(g(vef_b))
    return si
    
        
def g(v):
    return ((v+1)/2)*(np.log2((v+1)/2)) - ((v-1)/2)*(safe_log2((v-1)/2))    
    
    
def safe_log2(x, minval=0.0000000001):
    return np.log2(x.clip(min=minval))

def symplectic_eigenvalues(cm):
    N = int(cm.shape[0]/2)
    omega = symplectic_form(N)
    eigvals = np.linalg.eigvals(1j*np.dot(omega, cm)).real
    eigvals = eigvals[eigvals > 0]
    return eigvals


def symplectic_form(N):
    w = np.array([[0, 1], [-1, 0]])

    # Use direct sum to calculate the symplectic form
    return la.block_diag(*[w]*N)


def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)


def RCI_simple(rho, pos_keep):
    rho_a = rho.ptrace(pos_keep)
    rci = qt.entropy_vn(rho_a, base=2) - qt.entropy_vn(rho, base=2)
    return rci


def CI(sys, pos_trace_out):
    # Get both covariance matrices
    sys.set_quadratures_basis()
    CM1 = sys.get_full_CM()
    
    sys.ptrace(pos_trace_out)
    sys.set_quadratures_basis()
    CM2 = sys.get_full_CM()
    
#    vef1 = np.max(symplectic_eigenvalues(CM1))
#    vef2 = np.max(symplectic_eigenvalues(CM2))
    vef1 = symplectic_eigenvalues(CM1)
    vef2 = symplectic_eigenvalues(CM2)
    
    ci = np.sum(g(vef2)) - np.sum(g(vef1))
#    ci = sum(g(vef2)) - sum(g(vef1))
    return ci