import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op


def Delta(r, t, g, T, nth):
    g = g * T
    res = np.exp(-2*r - t)*((1 + np.exp(t/2) * g)**2 + np.exp(4*r) * (1 - np.exp(t/2) * g)**2 + \
          2 * np.exp(2*r + t) * (1 + g**2 + 2 * Gamma(r,t, g, T, nth)))
    return res


def Gamma(r, t, g, T, nth):
    R = 1 - T
    res = (1 - np.exp(-t)) * (.5 + nth) + g**2 * R**2
    return res


def squeezed_bell_eq(r, d, t, g, T, Br, Bi, nth):
    g = g * T
    Bnorm = Br**2 + Bi**2
    D = Delta(r, t, g, T, nth)
    res = (4/D) * np.exp((-4/D) * (g - 1)**2 * Bnorm) * (1 + (2*np.exp(-4*r - 2*t))/(D**4) * \
          ((1 + np.exp(t/2) * g ) - np.exp(4 * r) * (1 - np.exp(t/2) * g)**2)**2 * \
          (D**2 - 8 * D * (g -  1)**2 * Bnorm + 8 * (g - 1)**4 * Bnorm**2) * np.sin(d)**2 + \
          2 * np.exp(-2*r - 2*t)/D**2 * (4 * (g - 1)**2 * Bnorm - D) * np.sin(d) * \
          (np.cos(d) * (-(1 + np.exp(t/2)*g)**2 + np.exp(4*r) * (1 - np.exp(t/2)*g)**2) + \
          np.sin(d) * ((1 + np.exp(t/2)*g)**2 + np.exp(4*r) * (1 - np.exp(t/2)*g)**2)))
    return res


def squeezed_cat_eq2(r, d, gm, t, g, T, Br, Bi, nth):
    g = g * T
    Bnorm = Br**2 + Bi**2
    D = Delta(r, t, g, T, nth)
    # print('D', D)
    res = (4/(D * (1 + np.exp(-gm**2) * np.sin(2*d)))) * \
          (np.cos(d)**2 * np.exp(-(4/D) * (g-1)**2 * Bnorm) + \
          np.exp(-gm**2 - (1/D) * ((g-1) * 2 * Br - np.exp(r)*g*gm)**2) * \
          np.sin(d) * np.cos(d) * 2 * np.real(np.exp((1/D) * ((g-1)* 2j * Bi + np.exp(-r)*g*gm)**2)) + \
          np.sin(d)**2 * np.exp((-4/D) * (((g-1)*Br-np.exp(r)*gm*(g-np.exp(t/2)))**2 + ((g-1)*Bi)**2)))
    return res


def squeezed_cat_eq(r, d, gm, t, g, T, Br, Bi, nth):
    g = g * T
    Bnorm = Br**2 + Bi**2
    D = Delta(r, t, g, T, nth)
    pa = 4 / (D * (1 + np.exp(-gm**2) * np.sin(2*d)))
    pb = (-4 / D) * (g - 1)**2 * Bnorm
    pc = -gm**2 - (1 / D) * ((g - 1) * 2 * Br - np.exp(r) * g * gm)**2
    pd = 2 * np.real(np.exp((1 / D) * ((g - 1) * 2j * Bi + np.exp(-r) * g * gm)**2))
    pe = (-4 / D) * (((g-1) * Br - np.exp(r) * gm * (g - np.exp(t/2)))**2 + ((g-1) * Bi)**2)
    res = pa * (np.cos(d)**2 * np.exp(pb) + np.exp(pc) *np.sin(d) * np.cos(d) * pd + np.sin(d)**2 * np.exp(pe))
    print('res:', res)
    return res


def squeezed_cat_eq3( r, d, gm, t, g, T, Br, Bi, nth):
    g = g * T
    Bnorm = Br**2 + Bi**2
    D = Delta(r, t, g, T, nth)
    pa = 4 / (D * (1 + np.exp(-gm**2) * np.sin(2*d)))
    pb = (-4 / D) * (g - 1)**2 * Bnorm
    pc = -gm**2 - (1 / D) * ((g - 1) * 2 * Br - np.exp(r) * g * gm)**2
    pd = 2 * np.real(np.exp((1 / D) * ((g - 1) * 2j * Bi + np.exp(-r) * g * gm)**2))
    pe = (-4 / D) * (((g-1) * Br - np.exp(r) * gm * (g - np.exp(t/2)))**2 + ((g-1) * Bi)**2)
    res = pa * (np.cos(d)**2 * np.exp(pb) + np.exp(pc) *np.sin(d) * np.cos(d) * pd + np.sin(d)**2 * np.exp(pe))
    return res


def squeezed_cat_red(r, d, gm, nth):
    D = Delta(r, 0, 1, 1, nth)
    # print('D', D)
    a = 4 / (D * (1 + np.exp(-gm**2) * np.sin(2*d)))
    b = -gm**2 + gm**2/D * (np.exp(-2*r) - np.exp(2*r))
    res = a * (1 + 2 * np.cos(d) * np.sin(d) * np.exp(b))
    return res

def tmsv_eq(r, t, g, T, Br, Bi, nth):
    g = g * T
    Bnorm = Br**2 + Bi**2
    D = Delta(r, t, g, T, nth)
    # print('Delta:', D)
    res = (4/D) * np.exp((-4/D) * (g-1)**2 * Bnorm)
    return res


def pars_squeezed_bell(pars, t, g, T, Br, Bi, nth):
    r, d = pars
    return squeezed_bell_eq(r, d, t, g, T, Br, Bi, nth)

def pars_squeezed_cat(pars, t, g, T, Br, Bi, nth):
    r, d, gm = pars
    return squeezed_cat_eq(r, d, gm, t, g, T, Br, Bi, nth)

def pars_squeezed_bell_r(pars, r, t, g, T, Br, Bi, nth):
    d = pars[0]
    return squeezed_bell_eq(r, d, t, g, T, Br, Bi, nth)

def pars_squeezed_cat_r(pars, r, t, g, T, Br, Bi, nth):
    d, gm = pars
    # print(r, d, gm)
    # print(squeezed_cat_eq(r, d, gm, t, g, T, Br, Bi))
    return squeezed_cat_eq(r, d, gm, t, g, T, Br, Bi, nth)

def pars_tmsv(pars, t, g, T, Br, Bi, nth):
    return tmsv_eq(pars, t, g, T, Br, Bi, nth)


def opt_squeezed_bell(t, g, T, B, nth):
    Br = np.real(B)
    Bi = np.imag(B)
    F = lambda pars : 1 - pars_squeezed_bell(pars, t, g, T, Br, Bi, nth)
    initial_guess = [1, 1]
    res = op.minimize(F, initial_guess)
    print(res)
    return res['x'], res['success']


def opt_squeezed_cat(t, g, T, B, nth):
    Br = np.real(B)
    Bi = np.imag(B)
    F = lambda pars : 1 - pars_squeezed_cat(pars, t, g, T, Br, Bi, nth)
    initial_guess = [1, 1, 1]
    res = op.minimize(F, initial_guess)
    return res['x'], res['success']


def opt_tmsv(t, g, T, B, nth):
    Br = np.real(B)
    Bi = np.imag(B)
    F = lambda pars : 1 - pars_tmsv(pars, t, g, T, Br, Bi, nth)
    initial_guess = 1
    cons=({'type': 'ineq',
           'fun': lambda x: x})
    res = op.minimize(F, initial_guess)
    print(res)
    return res['x'], res['success']


def opt_squeezed_bell_r(r, t, g, T, B, nth):
    Br = np.real(B)
    Bi = np.imag(B)
    F = lambda pars : 1 - pars_squeezed_bell_r(pars, r, t, g, T, Br, Bi, nth)
    initial_guess = [1]
    res = op.minimize(F, initial_guess)
    print(res)
    return res['x'], res['success']


def opt_squeezed_cat_r(r, t, g, T, B, nth):
    Br = np.real(B)
    Bi = np.imag(B)
    F = lambda pars : 1 - pars_squeezed_cat_r(pars, r, t, g, T, Br, Bi, nth)
    initial_guess = [1, 1]
    res = op.minimize(F, initial_guess)
    return res['x'], res['success']


def get_F_func(stype):
    if stype == 'squeezed_bell':
        return pars_squeezed_bell
    elif stype == 'squeezed_cat':
        return pars_squeezed_cat
    elif stype == 'tmsv':
        return pars_tmsv

def get_F_func_r(stype):
    if stype == 'squeezed_bell':
        return pars_squeezed_bell_r
    elif stype == 'squeezed_cat':
        return pars_squeezed_cat_r
    elif stype == 'tmsv':
        return pars_tmsv_r

def get_opt_func(stype):
    if stype == 'squeezed_bell':
        return opt_squeezed_bell
    elif stype == 'squeezed_cat':
        return opt_squeezed_cat
    elif stype == 'tmsv':
        return opt_tmsv


def get_opt_func_r(stype):
    if stype == 'squeezed_bell':
        return opt_squeezed_bell_r
    elif stype == 'squeezed_cat':
        return opt_squeezed_cat_r


def get_opt_f(stype, t, g, T, B, nth):
    F_func = get_F_func(stype)
    opt_func = get_opt_func(stype)
    sol, success = opt_func(t, g, T, B, nth)
    if not success:
        print(sol)
        # raise AssertionError('Failure in optimization')
    return F_func(sol, t, g, T, np.real(B), np.imag(B), nth)


def get_opt_f_r(stype, r, t, g, T, B, nth):
    if stype == 'tmsv':
        return tmsv_eq(r, t, g, T, np.real(B), np.imag(B), nth)
    F_func = get_F_func_r(stype)
    opt_func = get_opt_func_r(stype)
    sol, success = opt_func(r, t, g, T, B, nth)
    if not success:
        print(sol)
        # raise AssertionError('Failure in optimization r')
    return F_func(sol, r, t, g, T, np.real(B), np.imag(B), nth)
