# -*- coding: utf-8 -*-
"""
Basic CV operations

Created on Fri Feb 15 12:05:43 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import tools
import hamiltonians



def displace(a, alpha):
    return (alpha*a.dag() - alpha.conjugate()*a).expm()


def squeeze(N, r, pos=0, Nmodes=1):
#    return ((z.conjugate()*a**2 - z*(a.dag()**2))/2).expm()
    # TODO: check if phase of pi is required
    S = qt.squeeze(N, -r)
    S = tools.tensor(N, S, pos, Nmodes)
    return S


def tmsqueeze(N, r, pos=[0,1], Nmodes=2):
    a = qt.tensor(qt.destroy(N), qt.identity(N))
    b = qt.tensor(qt.identity(N), qt.destroy(N))
    
#    S = (z.conjugate()*a*b - z*a.dag()*b.dag()).expm()
    # TODO: check this factor of two, and phase of pi
    S = qt.squeezing(a, b, -2*r)
    if Nmodes > 2:
        S = tools.reorder_two_mode_operator(N, S, pos, Nmodes)
    return S    


def beam_splitter(N, z, pos=[0,1], Nmodes=2):
    print("BS:", z)
    H = hamiltonians.beam_spliter(N, z)
    U = (-1j * H).expm()
    
    if Nmodes > 2:
        U = tools.reorder_two_mode_operator(N, U, pos, Nmodes)
    return U


def homodyne_operator2(N, phase, amplitude=1):
    z = amplitude*np.exp(1j*phase)
    print(z)
#    S = 1j*z*qt.create(N) - 1j*z.conjugate()*qt.destroy(N)
#    S = qt.squeeze(N, z)
    
    X = (qt.destroy(N) + qt.create(N))
    Y = -1j*(qt.destroy(N) - qt.create(N))
    S = abs(z)*(X*np.sin(np.angle(z)) + Y*np.cos(np.angle(z)))
    
    return S


def homodyne_operator(a, phase, amplitude=1):
    z = amplitude*np.exp(1j*phase)
#    S = 1j*z*qt.create(N) - 1j*z.conjugate()*qt.destroy(N)
#    S = qt.squeeze(N, z)
    
    X = (a + a.dag())
    Y = -1j*(a - a.dag())
    S = abs(z)*(X*np.sin(np.angle(z)) + Y*np.cos(np.angle(z)))
    
    return S


def var_homodyne(state, phase, amplitude=1):
    Shd = homodyne_operator2(state.dims[0][0], phase, amplitude)
    
    return qt.expect(Shd*Shd, state) - qt.expect(Shd, state)**2
    

def mean_homodyne(state, phase, amplitude=1):
    Shd = homodyne_operator2(state.dims[0][0], phase, amplitude)
    
    return qt.expect(Shd, state)


def purify(rho):
    eigenvals, eigenstates = rho.eigenstates()
    
    for i in range(len(eigenvals)):
        eigenstates[i] = np.sqrt(eigenvals[i]) * qt.tensor(eigenstates[i], eigenstates[i])
        
    return sum(eigenstates)


def photon_on_projector(N):
    P = 0
    for i in range(1, N):
        P += qt.basis(N, i).dag()
    return P 


#def photon_number_projector(n, N):
#    P = qt.basis(N, n).dag()
#    return P


def p_detect_photon_number_dm(N, rho, n, pos, Nmodes):
    P = qt.basis(N, n) * qt.basis(N, n).dag()
    P = tools.tensor(N, P, pos, Nmodes)
    return (P * rho).tr()


def p_detect_photon_number_ket(N, rho, n, pos, Nmodes):
    P = qt.basis(N, n).dag()
    P = tools.tensor(N, P, pos, Nmodes)
    return (P * rho).norm()


def collapse_photon_number(N, rho, n, pos, Nmodes):
    if rho.isket:
        return collapse_photon_number_ket(N, rho, n, pos, Nmodes)


def  collapse_photon_number_ket(N, rho, n, pos, Nmodes):
    P = qt.basis(N, n).dag()
    P = tools.tensor(N, P, pos, Nmodes)
    p = p_detect_photon_number_ket(N, rho, n, pos, Nmodes)
    rho = (P * rho)/p
    return rho, p


def  collapse_photon_number_dm(N, rho, n, pos, Nmodes):
    P = qt.basis(N, n).dag()
    P = tools.tensor(N, P, pos, Nmodes)
    rho = (P * rho)/p_detect_photon_number_dm(N, rho, n, pos, Nmodes)
    return rho


def RCI(rho, pos_keep):
    rho_a = rho.ptrace(pos_keep)
    rci = qt.entropy_vn(rho_a) - qt.entropy_vn(rho)
    return rci
    



