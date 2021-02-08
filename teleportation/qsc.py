# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:18:53 2021

@author: z5239621
"""

import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt

def fidelity(V, T, tsc, eps, eta, g, alpha):
    print('pars opt:', np.round(V, 3), np.round(tsc, 3), np.round(g, 3))

    gsc = np.sqrt((1 - tsc)/tsc)
    g = g*eta
    gt = g - 1

#    print(gsc, gt)

    A = 1/2 * (T*(V-1) + eps + 2)
    A2 = 1/ ((1 + gsc**2) *A)
    A3 = T * (V**2 - 1)/ (4*A)
    B1 = (gsc / (2 * A)) * np.sqrt(T * (V**2 - 1))
    A4 = 1/2 * (V + 2*A3 + g**2)
    A5 = 1 - gsc**2 * (1 - 1/A)
    A6 = -B1*g + gsc**2 * (A3/A - g * (1 - 1/A))
    A7 = (gsc * g)**2 * A3/A
    A8 = g**2/2 + A4 + g**2*(1-eta**2)/2
    A9 = A8 + 1/2
    A10 = gt**2/A9


    # print('A:', A)
    # print('A2:', A2)
    # print('A3:', A3)
    # print('B1:', B1)
    # print('A4:', A4)
    # print('A5:', A5)
    # print('A6:', A6)
    # print('A7:', A7)
    # print('A8:', A8)
    # print('A9:', A9)
    # print('A10:', A10)

    F = A2*np.exp(-A10 * np.abs(alpha)**2) * (A5/A9 + ((A6/(np.sqrt(2)*A9**3)) \
        * (4*A9 - 4*gt**2*np.abs(alpha)**2)) + ((A7/(2*A9**5)) \
        * (28*A9**2 - 56*A9*(gt*np.abs(alpha))**2 + 16*gt**4*(np.real(alpha)**4) \
           + np.imag(alpha)**4 + (np.real(alpha)*np.imag(alpha))**2)))



    F2 = A2 * np.exp(-A10*np.abs(alpha)**2) * (A5/A9 + (A6*2*np.sqrt(2))/(A9**2) + 14*A7/(A9**3) \
                      - (np.abs(alpha)**2 * ((2*np.sqrt(2)*g**2*A6)/A9**3 + (28*g**2*A7)/A9**4)) \
                      + 8*g**4*A7/A9**5 * (np.real(alpha)*np.imag(alpha))**2 \
                      + 8*g**4*A7/A9**5 * (np.real(alpha)**4 + np.imag(alpha)**4))

    return F


def fidelity_pars(pars, T, eps, eta, alpha):
    V, tsc, g = pars
    return fidelity(V, T, tsc, eps, eta, g, alpha)

def opt_fidelity(T, eps, eta, alpha):
    F = lambda P : 1 - fidelity_pars(P, T, eps, eta, alpha)
    initial_guess = [2, .2, 1]
    cons=({'type': 'ineq',
           'fun': lambda x: x[0] - 1},
          {'type': 'ineq',
           'fun': lambda x: x[1]},
          {'type': 'ineq',
           'fun': lambda x: 1.1-x[1]},
          {'type': 'ineq',
           'fun': lambda x: x[2]})
    

    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
    print(res)
#    print('opt V:', np.round(res['x'],3))
    return fidelity_pars(res['x'], T, eps, eta, alpha)


################################################
V = np.random.rand()*4 + 1
T = np.random.rand()
tsc = np.random.rand() * 0.5
tsc = 0.2
T = .005
V = 1.8
eps = 0
alpha = 1

print('pars:', V, T, tsc)
print('F:', fidelity(V, T, tsc, eps, 1, 1, alpha))

print('F_opt', opt_fidelity(T, eps, 1, alpha))

#print('F hand:', fidelity_pars([1.6, .999, .541], T, eps, 1, alpha))