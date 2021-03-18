# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:16:12 2021

@author: z5239621
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op


def squeezed_bell_avg(V, T, d, eps, eta, g, sigma):
    tau = -np.log(T)
    r = np.arccosh(V)/2
    nth = eps/((1-T)*2)
    return squeezed_bell_eq_avg(r, d, tau, g, eta, sigma, nth)


def squeezed_bell(V, T, d, eps, eta, g, alpha):
    tau = -np.log(T)
    r = np.arccosh(V)/2
    nth = eps/((1-T)*2)
    return squeezed_bell_eq(r, d, tau, g, eta, np.real(alpha), np.imag(alpha), nth)

def squeezed_bell_avg_pars(P, T, eps, eta, sigma):
    V, d, g = P
    return squeezed_bell_avg(V, T, d, eps, eta, g, sigma)


def opt_avg_fidelity(T, eps, eta, sigma):
    F = lambda P : 1 - squeezed_bell_avg_pars(P, T, eps, eta, sigma)
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
    return squeezed_bell_avg_pars(res['x'], T, eps, eta, sigma)



def squeezed_bell_eq(r, d, t, g, T, Br, Bi, nth):
    g = g * T
    Bnorm = Br**2 + Bi**2
    D = Delta(r, t, g, T, nth)
    res = (4/D) * np.exp((-4/D) * (g - 1)**2 * Bnorm) * (1 + (2*np.exp(-4*r - 2*t))/(D**4) * \
          ((1 + np.exp(t/2) * g )**2 - np.exp(4 * r) * (1 - np.exp(t/2) * g)**2)**2 * \
          (D**2 - 8 * D * (g -  1)**2 * Bnorm + 8 * (g - 1)**4 * Bnorm**2) * np.sin(d)**2 + \
          2 * np.exp(-2*r - 2*t)/D**2 * (4 * (g - 1)**2 * Bnorm - D) * np.sin(d) * \
          (np.cos(d) * (-(1 + np.exp(t/2)*g)**2 + np.exp(4*r) * (1 - np.exp(t/2)*g)**2) + \
          np.sin(d) * ((1 + np.exp(t/2)*g)**2 + np.exp(4*r) * (1 - np.exp(t/2)*g)**2)))
    return res

def squeezed_bell_eq_avg(r, d, t, g, T, sigma, nth):
    g = g * T
    D = Delta(r, t, g, T, nth)
    
    A1 = np.sin(d)**2 * 2 *np.exp(-4*r - 2*t)/D**4 * ((1+np.exp(t/2)*g)**2 - \
                np.exp(4*r)*(1-np.exp(t/2)*g)**2)**2 
    A2 = 2 * np.exp(-2*r - t)/D**2 * np.sin(d) * (np.cos(d)*(-(1+np.exp(t/2)*g)**2 + \
                                            np.exp(4*r) * ((1-np.exp(t/2)*g))**2) + \
                np.sin(d)* ((1+np.exp(t/2)*g)**2 + np.exp(4*r) * (1-np.exp(t/2)*g)**2))
    
    A = 1 + A1*D**2 - A2*D
    B = - A1 * 8 * D * (g-1)**2 + A2 * 4 * (g-1)**2
    C = A1 * 8 * (g-1)**4
    E = 4/D * (g-1)**2 + 1/sigma
    
    res = 4/(D * sigma) * (A/E + B/E**2 + C/E**3 * (12 + 1/2))
    return res

def Delta(r, t, g, T, nth):
    g = g * T
    res = np.exp(-2*r - t)*((1 + np.exp(t/2) * g)**2 + np.exp(4*r) * (1 - np.exp(t/2) * g)**2 + \
          2 * np.exp(2*r + t) * (1 + g**2 + 2 * Gamma(r,t, g, T, nth)))
    return res

def Gamma(r, t, g, T, nth):
    R = 1 - T
    res = (1 - np.exp(-t)) * (.5 + nth) + g**2 * R**2
    return res