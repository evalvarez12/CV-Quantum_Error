# -*- coding: utf-8 -*-
"""
Basic CV operations

Created on Fri Feb 15 12:05:43 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np

def beam_splitter(operators, theta):
    """
    Returns the corresponding beam splitter operator.
    Asummes symmetric beam splitter
    
    Returns
    --------
    bs_ops : qobj - Operator corresponding to the beam splitter transform
    """
#    U = np.array([[np.cos(theta) + 1j*np.sin(theta)] , [1j*np.sin(theta) + np.cos(theta)]])
    
    a, b = operators
    
    N = a.dims[0][0]
    identity = qt.identity(N)
    a = qt.tensor(a, identity)
    b = qt.tensor(identity, b)
    
    operators_out = [a*np.cos(theta) + b*1j*np.sin(theta), a*1j*np.sin(theta) + b*np.cos(theta)]
    return operators_out


def beam_splitter_inv(operators, theta):
    """
    Returns the corresponding beam splitter operator.
    Asummes symmetric beam splitter
    
    Returns
    --------
    bs_ops : qobj - Operator corresponding to the beam splitter transform
    """
#    U = np.array([[np.cos(theta) + 1j*np.sin(theta)] , [1j*np.sin(theta) + np.cos(theta)]])
    
    a, b = operators
    operators_out = [a*np.cos(theta) - b*1j*np.sin(theta), -a*1j*np.sin(theta) + b*np.cos(theta)]
    return operators_out


def tritter(operators, theta1, theta2):
    a, b, c = operators
    a, b = beam_splitter([a, b], theta1)
    b, c = beam_splitter([b, c], theta2)
    
    return [a, b, c]


def loss_channel(eta, a_in):
    theta = np.arccos(np.sqrt(eta))
    
    a_noise =qt.destroy(a_in.dims[0][0])
    
    a_out, _ = beam_splitter([a_in, a_noise], theta)
    return a_out


def homodyne_operator(N, phase, amplitude=1):
    z = amplitude*np.exp(1j*phase)
    print(z)
#    S = 1j*z*qt.create(N) - 1j*z.conjugate()*qt.destroy(N)
#    S = qt.squeeze(N, z)
    
    X = (qt.destroy(N) + qt.create(N))
    Y = -1j*(qt.destroy(N) - qt.create(N))
    S = abs(z)*(X*np.sin(np.angle(z)) + Y*np.cos(np.angle(z)))
    
    return S


def S_homodyne(state, phase, amplitude=1):
    Shd = homodyne_operator(state.dims[0][0], phase, amplitude)
    
    return qt.expect(Shd*Shd, state) - qt.expect(Shd, state)**2
    