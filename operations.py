# -*- coding: utf-8 -*-
"""
Basic CV operations

Created on Fri Feb 15 12:05:43 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import tools



def displace(a, alpha):
    return (alpha*a.dag() - alpha.conjugate()*a).expm()

def squeeze(a, z):
    return ((z.conjugate()*a**2 - z*(a.dag()**2))/2).expm()

def beam_splitter(operators, theta, tensor=True):
    """
    Returns the corresponding beam splitter operator.
    Asummes symmetric beam splitter
    
    Returns
    --------
    bs_ops : qobj - Operator corresponding to the beam splitter transform
    """
#    U = np.array([[np.cos(theta) + 1j*np.sin(theta)] , [1j*np.sin(theta) + np.cos(theta)]])
    
    a, b = operators
    
    if tensor:
        # Tensoring additional identities required to form the new Hlbert space
        N = a.dims[0][0]
        identity = qt.identity(N)
        a = qt.tensor(a, identity)
        b = qt.tensor(identity, b)
    
    operators_out = [a*np.cos(theta) + b*1j*np.sin(theta), a*1j*np.sin(theta) + b*np.cos(theta)]
    return operators_out


def tritter(operators, theta1, theta2, tensor=True):
    """
           |    
    c ->---/----->
           |
           |
    b ->---/----->
           |
           ^
           a
    """
    a, b, c = operators
    a, b = beam_splitter([a, b], theta1, tensor)
    
    if tensor:
        # Tensor additional identityes reqiured
        c = qt.tensor(qt.identity(c.dims[0][0]), c)
    
    a, c = beam_splitter([a, c], theta2, tensor)
    
    if tensor:
        b = qt.tensor(b, qt.identity(b.dims[0][0]))
    
    return [a, b, c]


def beam_splitter_Uoperator(N, theta, pos=[0,1]):
    Nextra = max(pos)
#    a = qt.tensor(qt.destroy(N), qt.identity(N))
    a = tools.tensor(qt.destroy(N), N, pos[0], Nextra)
    b = tools.tensor(qt.destroy(N), N, pos[1], Nextra)
#    b = qt.tensor(qt.identity(N), qt.destroy(N))
    
    
    # Define the mode-mixing Hamiltonian
    H = (a.dag()*b + a*b.dag())
    
    # Unitary evolution
    U = (-1j*theta*H).expm()
    
    return U


def beam_splitter_applyU(rho, theta):
    # Get unitary opertator
    U = beam_splitter_Uoperator(rho.dims[0][0], theta)
    
    rho = U * rho * U.dag()
    return rho


def tritter_applyU(rho, theta1, theta2):
    """
           |    
    c ->---/----->
           |
           |
    b ->---/----->
           |
           ^
           a
    """
    N = rho.dims[0][0]
    
    U1 = beam_splitter_Uoperator(N, theta1)
    U1 = qt.tensor(U1, qt.identity(N))
    U2 = beam_splitter_Uoperator(N, theta2)
    U2 = qt.tensor(qt.identity(N), U2)

    rho = U1 * rho * U1.dag()
    rho = U2 * rho * U2.dag()




def loss_channel(eta, a_in):
    theta = np.arccos(np.sqrt(eta))
    
    a_noise =qt.destroy(a_in.dims[0][0])
    
    a_out, _ = beam_splitter([a_in, a_noise], theta)
    return a_out


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

