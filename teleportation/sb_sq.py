# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:16:12 2021

@author: z5239621
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import scipy.integrate as integrate

def squeezed_bell_eq(r, d, t, g, T, s, nth):
    G = Gamma(r, t, g, T, nth)
    g2 = g * T
    D1 = 1 + np.exp(4*r) + 2*np.exp(t/2)*(1-np.exp(4*r))*g2 + np.exp(t)*(1+np.exp(4*r))*g2**2
    D2 = 1 - np.exp(4*r) + 2*np.exp(t/2)*(1+np.exp(4*r))*g2 + np.exp(t)*(1-np.exp(4*r))*g2**2
    L1 = np.exp(-2*r - t)*D1 + 2*np.exp(2*s)*(1+g2**2) + 4*G
    L2 = np.exp(-2*r - t)*D1 + 2*np.exp(-2*s)*(1+g2**2) + 4*G

    # pp = np.array([np.round(L1,3), np.round(L2, 3), np.round(G, 3), np.round(D1, 3), np.round(D2, 3)])
    # if not np.isnan(pp).any():
    #     print(pp)

    si = np.sin(d)
    co = np.cos(d)

    res = 4/np.sqrt(L1*L2)*(1 + np.exp(-2*r-t)*si*(D2*co-D1*si)*((1/L1) + (1/L2)) + \
                    1/4*np.exp(-2*(2*r+t))*(D2*si)**2*((3/L1**2) + (3/L2**2) + \
                    2/(L1*L2)))

    return res


def squeezed_bell_eq_disp(r, d, t, g, T, a, s, nth):
    G = Gamma(r, t, g, T, nth)
    g2 = g * T
    D1 = 1 + np.exp(4*r) + 2*np.exp(t/2)*(1-np.exp(4*r))*g2 + np.exp(t)*(1+np.exp(4*r))*g2**2
    D2 = 1 - np.exp(4*r) + 2*np.exp(t/2)*(1+np.exp(4*r))*g2 + np.exp(t)*(1-np.exp(4*r))*g2**2
    L1 = np.exp(-2*r - t)*D1 + 2*np.exp(2*s)*(1+g2**2) + 4*G
    L2 = np.exp(-2*r - t)*D1 + 2*np.exp(-2*s)*(1+g2**2) + 4*G
    w1 = (1-g2)**2*(2*np.imag(a))**2
    w2 = (1-g2)**2*(2*np.real(a))**2
    # pp = np.array([np.round(L1,3), np.round(L2, 3), np.round(G, 3), np.round(D1, 3), np.round(D2, 3)])
    # if not np.isnan(pp).any():
    #     print(pp)

    si = np.sin(d)
    co = np.cos(d)

    E = np.exp(w1**2/L1 - w2**2/L2)
    res = 4/np.sqrt(L1*L2)*E*(1 + np.exp(-2*r-t)*si*(D2*co-D1*si)*((1/L1)*(1 + 2*w1**2/L1) + (1/L2)*(1 + 2*w2**2/L2)) + \
                    1/4*np.exp(-2*(2*r+t))*(D2*si)**2*((1/L1**2)*(3 + 12*w1**2/L1 + 4*w1**4/L1**2) + (1/L2**2)*(3 + 12*w2**2/L2 + 4*w2**4/L2**2) + \
                    (2/(L1*L2))*(1 + 2*w1**2/L1 - 2*w2**2/L2 - 4*w1**2*w2**2/(L1*L2))))

    return res


def Gamma(r, t, g, T, nth):
    # NOTE: this fix of /2 to R
    R = (1 - T**2)/2
    res = (1 - np.exp(-t)) * (.5 + nth) + g**2 * R
    return res



def avg_fidelity(V, T, d, eps, eta, g, sig):
    F = lambda s : fidelity(V, T, d, eps, eta, g, s) * np.exp(-s**2/sig) * s
    I = integrate.quad(F, 0, 3*sig)
#    print(I)
    return 2*I[0]/sig


def fidelity(V, T, d, eps, eta, g, s):

    tau = -np.log(T)
    r = np.arccosh(V)/2
    if T == 1:
        nth = eps
    else:
        nth = eps/((1-T)*2)

#    P = np.array([V, T, d, eps, eta, g, s])
#    P2 = np.array([r, d, tau, g, eta, s, nth])
    # if not np.isnan(P).any():
        # print(P)
        # print(P2)

    return squeezed_bell_eq(r, d, tau, g, eta, s, nth)

def fidelity_disp(V, T, d, eps, eta, g, a, s):
    tau = -np.log(T)
    r = np.arccosh(V)/2
    if T == 1:
        nth = eps
    else:
        nth = eps/((1-T)*2)

#    P = np.array([V, T, d, eps, eta, g, s])
#    P2 = np.array([r, d, tau, g, eta, s, nth])
    # if not np.isnan(P).any():
        # print(P)
        # print(P2)

    return squeezed_bell_eq_disp(r, d, tau, g, eta, a, s, nth)


def fidelity_pars(P, T, eps, eta, s):
    V, d, g = P
    # if not np.isnan(P).any():
        # print(P)

    return fidelity(V, T, d, eps, eta, g, s)


def avg_fidelity_pars(P, T, eps, eta, sig):
    V, d, g = P
    # if not np.isnan(P).any():
        # print(P)

    return avg_fidelity(V, T, d, eps, eta, g, sig)


def fidelity_disp_pars(P, T, eps, eta, a, s):
    V, d, g = P
    # if not np.isnan(P).any():
        # print(P)

    return fidelity_disp(V, T, d, eps, eta, g, a, s)


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
           'fun': lambda x: x[2] - eta})

    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
#    print(res)
#    print('opt V:', np.round(res['x'],3))
    return fidelity_pars(res['x'], T, eps, eta, s)


def opt_avg_fidelity(T, eps, eta, sig):
    F = lambda P : 1 - avg_fidelity_pars(P, T, eps, eta, sig)
    initial_guess = [1.5, np.pi/3, 1]
    cons=({'type': 'ineq',
           'fun': lambda x: x[0] - 1.001},
          {'type': 'ineq',
           'fun': lambda x: x[1]},
          {'type': 'ineq',
           'fun': lambda x: 2*np.pi - x[1]},
          {'type': 'ineq',
           'fun': lambda x: x[2] - eta})

    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
#    print(res)
#    print('opt V:', np.round(res['x'],3))
    return avg_fidelity_pars(res['x'], T, eps, eta, sig)


def opt_fidelity_disp(T, eps, eta, a, s):
    F = lambda P : 1 - fidelity_disp_pars(P, T, eps, eta, a, s)
    initial_guess = [1.5, np.pi/3, 1.1]
    cons=({'type': 'ineq',
           'fun': lambda x: x[0] - 1.001},
          {'type': 'ineq',
           'fun': lambda x: 5 - x[0]},
          {'type': 'ineq',
           'fun': lambda x: x[1]},
          {'type': 'ineq',
           'fun': lambda x: 2*np.pi - x[1]},
          {'type': 'ineq',
           'fun': lambda x: x[2] - eta},
          {'type': 'ineq',
           'fun': lambda x: 2- x[2]})

    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
    # print(res)
#    print('opt V:', np.round(res['x'],3))
    return fidelity_disp_pars(res['x'], T, eps, eta, a, s)


def opt_fidelity_params(T, eps, eta, s):
    F = lambda P : 1 - fidelity_pars(P, T, eps, eta, s)
    initial_guess = [1.5, np.pi/3, 1.1]
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
    print(res)
#    print('opt V:', np.round(res['x'],3))
#    return fidelity_pars(res['x'], T, eps, eta, s)
    return res['x']


def opt_fidelity_avg(T, eps, eta, s_mean, sigma):
    s = np.random.normal(s_mean, np.sqrt(sigma), 20000)
    opt_params = opt_fidelity_params(T, eps, eta, s_mean)
    V_opt, d_opt, g_opt = opt_params

    F = fidelity(V_opt, T, d_opt, eps, eta, g_opt, s)
#    print(s)
#    print(V_opt, d_opt, g_opt)
#    print(F)
    return np.average(F)

def fidelity_pars_r(P, V, T, eps, eta, s):
    d, g = P

    return fidelity(V, T, d, eps, eta, g, s)

def opt_fidelity_r(V, T, eps, eta, s):
    F = lambda P : 1 - fidelity_pars_r(P, V, T, eps, eta, s)
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


# T = 0.99
# eta = .8
# eps = 0
# s = 5
#
# print(opt_fidelity(T, eps, eta, s))
