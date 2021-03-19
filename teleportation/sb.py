# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:16:12 2021

@author: z5239621
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op


def avg_fidelity(V, T, d, eps, eta, g, sigma):
    tau = -np.log(T)
    r = np.arccosh(V)/2
    nth = eps/((1-T)*2)
    return squeezed_bell_eq_avg(r, d, tau, g, eta, sigma, nth)


def fidelity(V, T, d, eps, eta, g, alpha):
    tau = -np.log(T)
    r = np.arccosh(V)/2
    nth = eps/((1-T)*2)
    return squeezed_bell_eq(r, d, tau, g, eta, np.real(alpha), np.imag(alpha), nth)


def avg_fidelity_pars_r(P, V, T, eps, eta, sigma):
    d, g = P
    return avg_fidelity(V, T, d, eps, eta, g, sigma)


def avg_fidelity_pars(P, T, eps, eta, sigma):
    V, d, g = P
    return avg_fidelity(V, T, d, eps, eta, g, sigma)


def fidelity_pars_r(P, V, T, eps, eta, alpha):
    d, g = P
    return fidelity(V, T, d, eps, eta, g, alpha)


def opt_avg_fidelity(T, eps, eta, sigma):
    F = lambda P : 1 - avg_fidelity_pars(P, T, eps, eta, sigma)
    initial_guess = [1.5, np.pi/3, 1]
    cons=({'type': 'ineq',
           'fun': lambda x: x[0] - 1.001},
          {'type': 'ineq',
           'fun': lambda x: 10000 - x[0]},
          {'type': 'ineq',
           'fun': lambda x: x[1]},
          {'type': 'ineq',
           'fun': lambda x: 2*np.pi - x[1]},
          {'type': 'ineq',
           'fun': lambda x: x[2]}) 

    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
#    print(res)
#    print('opt V:', np.round(res['x'],3))
    return avg_fidelity_pars(res['x'], T, eps, eta, sigma)


def opt_avg_fidelity_r(V, T, eps, eta, sigma):
    F = lambda P : 1 - avg_fidelity_pars_r(P, V, T, eps, eta, sigma)
    initial_guess = [np.pi/3, 1]
    cons=({'type': 'ineq',
           'fun': lambda x: x[0]},
          {'type': 'ineq',
           'fun': lambda x: 2*np.pi - x[0]},
          {'type': 'ineq',
           'fun': lambda x: x[1]}) 

    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
#    print(res)
#    print('opt V:', np.round(res['x'],3))
    return avg_fidelity_pars_r(res['x'], V, T, eps, eta, sigma)


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
    print(res)
#    print('opt V:', np.round(res['x'],3))
    return fidelity_pars_r(res['x'], V, T, eps, eta, alpha)


def squeezed_bell_eq(r, d, t, g, T, Br, Bi, nth):
    Bnorm = Br**2 + Bi**2
    D = Delta(r, t, g, T, nth)
    g = g * T
    
    res = (4/D) * np.exp((-4/D) * (g - 1)**2 * Bnorm) * (1 + (2*np.exp(-4*r - 2*t))/(D**4) * \
          ((1 + np.exp(t/2) * g )**2 - np.exp(4 * r) * (1 - np.exp(t/2) * g)**2)**2 * \
          (D**2 - 8 * D * (g -  1)**2 * Bnorm + 8 * (g - 1)**4 * Bnorm**2) * np.sin(d)**2 + \
          2 * np.exp(-2*r - 2*t)/D**2 * (4 * (g - 1)**2 * Bnorm - D) * np.sin(d) * \
          (np.cos(d) * (-(1 + np.exp(t/2)*g)**2 + np.exp(4*r) * (1 - np.exp(t/2)*g)**2) + \
          np.sin(d) * ((1 + np.exp(t/2)*g)**2 + np.exp(4*r) * (1 - np.exp(t/2)*g)**2)))
    return res

def squeezed_bell_eq_avg(r, d, t, g, T, sigma, nth):
    D = Delta(r, t, g, T, nth) 
    g = g * T
    
    A1 = np.sin(d)**2 * 2 *np.exp(-4*r - 2*t)/D**4 * ((1+np.exp(t/2)*g)**2 - \
                np.exp(4*r)*(1-np.exp(t/2)*g)**2)**2 
    A2 = 2 * np.exp(-2*r - t)/D**2 * np.sin(d) * (np.cos(d)*(-(1+np.exp(t/2)*g)**2 + \
                                            np.exp(4*r) * ((1-np.exp(t/2)*g))**2) + \
                np.sin(d)* ((1+np.exp(t/2)*g)**2 + np.exp(4*r) * (1-np.exp(t/2)*g)**2))
    
    A = 1 + A1*D**2 - A2*D
    B = - A1 * 8 * D * (g-1)**2 + A2 * 4 * (g-1)**2
    C = A1 * 8 * (g-1)**4
    D2 = 4/D * (g-1)**2 + 1/sigma
    
    res = 4/(D * sigma) * (A/D2 + B/D2**2 + C/D2**3 * (12 + 1/2))
    return res

def Delta(r, t, g, T, nth):
    G = Gamma(r, t, g, T, nth)
    g = g * T
    
    res = np.exp(-2*r - t)*((1 + np.exp(t/2) * g)**2 + np.exp(4*r) * (1 - np.exp(t/2) * g)**2 + \
          2 * np.exp(2*r + t) * (1 + g**2 + 2 * G))
    return res

def Gamma(r, t, g, T, nth):
    R = 1 - T
    res = (1 - np.exp(-t)) * (.5 + nth) + g**2 * R**2
    return res