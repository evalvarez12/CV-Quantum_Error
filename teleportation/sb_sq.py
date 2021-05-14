# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:16:12 2021

@author: z5239621
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op


def fidelity(V, T, d, eps, eta, g, s):
    tau = -np.log(T)
    r = np.arccosh(V)/2
    if T == 1:
        nth = eps
    else:
        nth = eps/((1-T)*2)
    return squeezed_bell_eq(r, d, tau, g, eta, s, nth)


def fidelity_pars_r(P, V, T, eps, eta, alpha):
    d, g = P
#    g = 1/eta
    return fidelity(V, T, d, eps, eta, g, alpha)


def fidelity_pars(P, T, eps, eta, s):
    V, d, g = P
#    g = 1/eta
    return fidelity(V, T, d, eps, eta, g, s)


def opt_fidelity_r(V, T, eps, eta, alpha):
    F = lambda P : 1 - fidelity_pars_r(P, V, T, eps, eta, alpha)
    initial_guess = [np.pi/3, 1]
    cons=({'type': 'ineq',
           'fun': lambda x: x[0]},
          {'type': 'ineq',
           'fun': lambda x: 2*np.pi - x[0]},
          {'type': 'ineq',
           'fun': lambda x: x[1]})

    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
    # print(res)
#    print('opt V:', np.round(res['x'],3))
    return fidelity_pars_r(res['x'], V, T, eps, eta, alpha)



def opt_fidelity(T, eps, eta, s):
    F = lambda P : 1 - fidelity_pars(P, T, eps, eta, s)
    initial_guess = [1.5, np.pi/3, 1]
    cons=({'type': 'ineq',
           'fun': lambda x: x[0] - 1.001},  
          {'type': 'ineq',
           'fun': lambda x: x[1]},
          {'type': 'ineq',
           'fun': lambda x: 2*np.pi - x[1]},
          {'type': 'ineq',
           'fun': lambda x: x[2]})

    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
    # print(res)
#    print('opt V:', np.round(res['x'],3))
    return fidelity_pars(res['x'], T, eps, eta, s)


def squeezed_bell_eq(r, d, t, g, T, s, nth):
    G = Gamma(r, t, g, T, nth)
    g2 = g * T
    D1 = 1 + np.exp(4*r) + 2*np.exp(t/2)*(1-np.exp(4*r))*g2 + np.exp(t)*(1+np.exp(4*r))*g2 
    D2 = 1 - np.exp(4*r) + 2*np.exp(t/2)*(1+np.exp(4*r))*g2 + np.exp(t)*(1-np.exp(4*r))*g2 
    L1 = np.exp(-2*r - t)*D1 + 2*np.exp(2*s)*(1+g2**2) + 4*G
    L2 = np.exp(-2*r - t)*D1 + 2*np.exp(-2*s)*(1+g2**2) + 4*G


    si = np.sin(d)
    co = np.cos(d)

    res = 4/np.sqrt(L1*L2)*(1 + np.exp(-2*r-t)*si*(D2*co-D1*si)*((1/L1) + (1/L2)) + \
                    1/4*np.exp(-2*(2*r+t))*(D2*si)**2*3*((1/L1**2) + (1/L2**2) + \
                    2/(L1*L2)))
        
    return res


def Gamma(r, t, g, T, nth):
    # NOTE: this fix of /2 to R
    R = (1 - T**2)/2 
    res = (1 - np.exp(-t)) * (.5 + nth) + g**2 * R
    return res
