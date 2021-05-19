# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:18:53 2021

@author: z5239621
"""

import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import TMSV as tmsv



def avg_fidelity(V, T, tsc, eps, eta, g, sigma):
    gsc = np.sqrt((1 - tsc)/tsc)
    g = g*eta
    gt = g - 1

#    print('gsc:', gsc)
#    print(np.round([V, tsc, g], 3))

    A = 0.5 * (T * (V - 1) + eps + 2)
    B = np.sqrt(T * (V**2 - 1))
    A2 = 0.5 * (V - (T * (V**2 - 1))/(T * (T - 1) + eps + 2))
    A3 = 1/A + gsc**2 * (1/A - 1/A**2)
    B1u = - gsc * B/(2*A**2) - B**2/4*A**3 - g**2 * (1/A - 1/A**2)
    B1v = gsc * B/(2*A**2) - B**2/4*A**3 - g**2 * (1/A - 1/A**2)
    A4 = g**2 * B**2 / (4 * A**3)
    A5 = A2 + 2 * g**2 + 1 - (g/eta)**2 * (1 - eta**2)
    A6 = A3/A5 + A4/A5**3 * (12 + 1/4) + 1/(2 * A5**2) * (B1u + B1v)
    B2u = - B1u * gt**2 / A5**3 - (A4 / A5**4) * gt**2 * (24 + 1/2)
    B2v = - B1v * gt**2 / A5**3 - (A4 / A5**4) * gt**2 * (24 + 1/2)
    A7 = (A4 / A5**5) * gt**4
    A8 = gt**2/A5 + 1/sigma
    
    norm = 1/A - gsc**2 / ((1 + gsc**2) * A**2)  
    
#    print('A1:', A1)
#    print('A2:', A2)
#    print('A3:', A3)
#    print('A5:', A5)
#    print('A6:', A6)
#    print('norm:', norm)

    F = (1 / (1 + gsc**2)) * 1/sigma * (A6/A8 + (1/(2*A8**2))* (B2u + B2v) + A7/(A8**3)* (96 + 1/4))
    return F
    

def fidelity(V, T, tsc, eps, eta, g, alpha):
#    print('pars opt:', np.round(V, 3), np.round(tsc, 3), np.round(g, 3))

    gsc = np.sqrt((1 - tsc)/tsc)
    g = g*eta
    gt = g - 1

    print('gsc:', gsc)
    print(np.round([V, tsc, g], 3))


    A = 0.5 * (T * (V - 1) + eps + 2)
    B = np.sqrt(T * (V**2 - 1))
    A2 = 0.5 * (V - (T * (V**2 - 1))/(T * (T - 1) + eps + 2))
    A3 = 1/A + gsc**2 * (1/A - 1/A**2)
    B1u = - gsc * B/(2*A**2) - B**2/4*A**3 - g**2 * (1/A - 1/A**2)
    B1v = gsc * B/(2*A**2) - B**2/4*A**3 - g**2 * (1/A - 1/A**2)
    A4 = g**2 * B**2 / (4 * A**3)
    A5 = A2 + 2 * g**2 + 1 - (g/eta)**2 * (1 - eta**2)
    A6 = A3/A5 + A4/A5**3 * (12 + 1/4) + 1/(2 * A5**2) * (B1u + B1v)
    B2u = - B1u * gt**2 / A5**3 - (A4 / A5**4) * gt**2 * (24 + 1/2)
    B2v = - B1v * gt**2 / A5**3 - (A4 / A5**4) * gt**2 * (24 + 1/2)
    A7 = (A4 / A5**5) * gt**4
    
    norm = 1/A - gsc**2 / ((1 + gsc**2) * A**2)

    
#    print('A1:', A1)
#    print('A2:', A2)
#    print('A3:', A3)
#    print('A4:', A4)
#    print('A5:', A5)
    print('norm:', norm)

    F = (1 / (1 + gsc**2)) * np.exp(-(gt**2/A5)*np.abs(alpha)**2) * \
        (A6 + B2u * np.imag(alpha)**2 + B2v * np.real(alpha)**2 + \
         A7 * (np.real(alpha)*np.imag(alpha))**2 + \
         8*A7 * (np.real(alpha)**4 + np.imag(alpha)**4))
    return F


def fidelity_pars(pars, T, eps, eta, alpha):
    V, tsc, g = pars
    return fidelity(V, T, tsc, eps, eta, g, alpha)


def avg_fidelity_pars(pars, T, eps, eta, sigma):
    V, tsc, g = pars
    return avg_fidelity(V, T, tsc, eps, eta, g, sigma)




def get_cons():
    cons=({'type': 'ineq',
       'fun': lambda x: x[0] - 1.001},
      {'type': 'ineq',
       'fun': lambda x: 10000 - x[0]},
      {'type': 'ineq',
       'fun': lambda x: x[1] - 0.001},
      {'type': 'ineq',
       'fun': lambda x: .99 - x[1]},
      {'type': 'ineq',
       'fun': lambda x: x[2]})
    return cons

def opt_fidelity(T, eps, eta, alpha):
    F = lambda P : 1 - fidelity_pars(P, T, eps, eta, alpha)
    initial_guess = [2, .2, 1]
    cons= get_cons()

    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
    print(res)
#    print('opt V:', np.round(res['x'],3))
    return fidelity_pars(res['x'], T, eps, eta, alpha)


def opt_avg_fidelity(T, eps, eta, sigma):
    F = lambda P : 1 - avg_fidelity_pars(P, T, eps, eta, sigma)
    initial_guess = [1.5, .2, 1]
    cons= get_cons()

    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
    print(res)
#    print('opt V:', np.round(res['x'],3))
    return avg_fidelity_pars(res['x'], T, eps, eta, sigma)

################################################
V = np.random.rand() + 1
T = np.random.rand() * 1
tsc = np.random.rand() 
#T = 1
#V = 1.00001
eps = 0.05
alpha = 5
sigma = 5
eta = 1
g = 1
#print('-------------Random')
#print('pars:', np.round([V, T, tsc], 3))
#print('F:', fidelity(V, T, tsc, eps, eta, g, alpha))
#print('F_avg:', avg_fidelity(V, T, tsc, eps, eta, g, sigma))
#print('-------------------')
#print('F_opt', opt_fidelity(T, eps, 1, alpha))

eps = 0.05
alpha = 5

sigma = 2
eta = 1
g = 1
T = 0.01
tsc = 0.04
V = 1.4
alpha = 10
gsc = np.sqrt((1 - tsc)/tsc)
#print('gsc:', gsc)
#print('p_succ:', 2 / (V+1) - (gsc*2/(V+1))**2 * (1 / (1 + gsc**2)))
#
#print('F:', fidelity(V, T, tsc, eps, eta, g, alpha))
#print('F_opt', opt_fidelity(T, eps, eta, sigma))
#print('F_TMSV:', tmsv.fidelity(V, T, eps, eta, g, alpha))
#print('F_avg:', avg_fidelity(V, T, .01, eps, eta, g, sigma))
#
#print('F_avg_TMSV:', tmsv.opt_fidelity_alphabet_gopt(T, eps, eta, sigma))
#print('F_avg_opt', opt_avg_fidelity(T, eps, eta, sigma))
#print('F hand:', fidelity_pars([1.6, .999, .541], T, eps, 1, alpha))