import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op


class Fidelity:
    def __init__(self):
        return
    # def __init__(self, rtype, beta, squeezing, tau, gain, T):
    #     self.rtype = rtype
    #     self.r = squeezing
    #     self.t = tau
    #     self.g = gain
    #     self.T = T
    #     self.B = beta
    #
    #     self.Delta = Delta(self.r, self.t, self.g)

    def Delta(self, r, t, g, T):
        G = self.Gamma(r, t, g, T, nth=0)
        g = g * T
        res = np.exp(-2*r - t)*((1 + np.exp(t/2) * g)**2 + np.exp(4*r) * (1 - np.exp(t/2) * g)**2 + \
              2 * np.exp(2*r + t) * (1 + g**2 + 2 * self.Gamma(r,t, g, T)))
        return res


    def Gamma(self, r, t, g, T, nth=0):
        R = 1 - T
        res = (1 - np.exp(-t)) * (.5 + nth) + g**2 * R**2
        return res


    def squeezed_bell_eq(self, r, d, t, g, T, Br, Bi):
        g = g * T
        Bnorm = Br**2 + Bi**2
        D = self.Delta(r, t, g, T)
        res = (4/D) * np.exp((-4/D) * (g - 1)**2 * Bnorm) * (1 + (2*np.exp(-4*r - 2*t))/(D**4) * \
              ((1 + np.exp(t/2) * g ) - np.exp(4 * r) * (1 - np.exp(t/2) * g)**2)**2 * \
              (D**2 - 8 * D * (g -  1)**2 * Bnorm + 8 * (g - 1)**4 * Bnorm**2) * np.sin(d)**2 + \
              2 * np.exp(-2*r - 2*t)/D**2 * (4 * (g - 1)**2 * Bnorm - D) * np.sin(d) * \
              (np.cos(d) * (-(1 + np.exp(t/2)*g)**2 + np.exp(4*r) * (1 - np.exp(t/2)*g)**2) + \
              np.sin(d) * ((1 + np.exp(t/2)*g)**2 + np.exp(4*r) * (1 - np.exp(t/2)*g)**2)))
        return res


    def squeezed_cat_eq(self, r, d, gm, t, g, T, Br, Bi):
        g = g * T
        Bnorm = Br**2 + Bi**2
        D = self.Delta(r, t, g, T)
        res = (4/(D * (1 + np.exp(-gm) * np.sin(2*d)))) * (np.cos(d)**2 * np.exp(-(4/D) * (g-1)**2 * Bnorm) + \
              np.exp(-gm - (1/D) * ((g-1) * 2 * Br - np.exp(r)*g*gm)**2) * \
              np.sin(d) * np.cos(d) * 2 * np.real(np.exp((1/D) * ((g-1)* 2j * Bi + np.exp(-r)*g*gm)**2)) + \
              np.sin(d)**2 * np.exp((-4/D) * (((g-1)*Br-np.exp(r)*gm*(g-np.exp(t/2)))**2 + Bi**2)))
        return res


    def tmsv_eq(self, r, t, g, T, Br, Bi):
        g = g * T
        Bnorm = Br**2 + Bi**2
        D = self.Delta(r, t, g, T)
        res = (4/D) * np.exp((-4/D) * (g-1)**2 * Bnorm)
        return res


    def pars_squeezed_bell(self, pars, t, g, T, Br, Bi):
        r, d = pars
        return self.squeezed_bell_eq(r, d, t, g, T, Br, Bi)

    def pars_squeezed_cat(self, pars, t, g, T, Br, Bi):
        r, d, gm = pars
        return self.squeezed_cat_eq(r, d, gm, t, g, T, Br, Bi)

    def pars_squeezed_bell_r(self, pars, r, t, g, T, Br, Bi):
        d = pars
        return self.squeezed_bell_eq(r, d, t, g, T, Br, Bi)

    def pars_squeezed_cat_r(self, pars, r, t, g, T, Br, Bi):
        d, gm = pars
        return self.squeezed_cat_eq(r, d, gm, t, g, T, Br, Bi)

    def pars_tmsv(self, pars, t, g, T, Br, Bi):
        return self.tmsv_eq(pars, t, g, T, Br, Bi)


    def opt_squeezed_bell(self, t, g, T, B):
        Br = np.real(B)
        Bi = np.imag(B)
        F = lambda pars : 1 - self.pars_squeezed_bell(pars, t, g, T, Br, Bi)
        initial_guess = [1, 1]
        res = op.minimize(F, initial_guess)
        return res['x'], res['success']


    def opt_squeezed_cat(self, t, g, T, B):
        Br = np.real(B)
        Bi = np.imag(B)
        F = lambda pars : 1 - self.pars_squeezed_cat(pars, t, g, T, Br, Bi)
        initial_guess = [1, 1, 1]
        res = op.minimize(F, initial_guess)
        return res['x'], res['success']


    def opt_tmsv(self, t, g, T, B):
        Br = np.real(B)
        Bi = np.imag(B)
        F = lambda pars : 1 - self.pars_tmsv(pars, t, g, T, Br, Bi)
        initial_guess = 1
        res = op.minimize(F, 1)
        return res['x'], res['success']


    def opt_squeezed_bell_r(self, r, t, g, T, B):
        Br = np.real(B)
        Bi = np.imag(B)
        F = lambda pars : 1 - self.pars_squeezed_bell_r(pars, r, t, g, T, Br, Bi)
        initial_guess = [1]
        res = op.minimize(F, initial_guess)
        return res['x'], res['success']


    def opt_squeezed_cat_r(self, r, t, g, T, B):
        Br = np.real(B)
        Bi = np.imag(B)
        F = lambda pars : 1 - self.pars_squeezed_cat_r(pars, r, t, g, T, Br, Bi)
        initial_guess = [1, 1]
        res = op.minimize(F, initial_guess)
        return res['x'], res['success']


    def get_F_func(self, stype):
        if stype == 'squeezed_bell':
            return self.pars_squeezed_bell
        elif stype == 'squeezed_cat':
            return self.pars_squeezed_cat
        elif stype == 'tmsv':
            return self.pars_tmsv

    def get_F_func_r(self, stype):
        if stype == 'squeezed_bell':
            return self.pars_squeezed_bell_r
        elif stype == 'squeezed_cat':
            return self.pars_squeezed_cat_r
        elif stype == 'tmsv':
            return self.pars_tmsv_r

    def get_opt_func(self, stype):
        if stype == 'squeezed_bell':
            return self.opt_squeezed_bell
        elif stype == 'squeezed_cat':
            return self.opt_squeezed_cat
        elif stype == 'tmsv':
            return self.opt_tmsv


    def get_opt_func_r(self, stype):
        if stype == 'squeezed_bell':
            return self.opt_squeezed_bell_r
        elif stype == 'squeezed_cat':
            return self.opt_squeezed_cat_r


    def get_opt_f(self, stype, t, g, T, B):
        F_func = self.get_F_func(stype)
        opt_func = self.get_opt_f(stype)
        sol, success = opt_func(t, g, T, B)
        if not success:
            raise AssertionError('Failure in optimization')
        return F_func(sol, t, g, T, np.real(B), np.imag(B))


    def get_opt_f_r(self, stype, r, t, g, T, B):
        if stype == 'tmsv':
            return self.tmsv_eq(r, t, g, T, np.real(B), np.imag(B))
        F_func = self.get_F_func_r(stype)
        opt_func = self.get_opt_func_r(stype)
        sol, success = opt_func(r, t, g, T, B)
        if not success:
            raise AssertionError('Failure in optimization')
        return F_func(sol, r, t, g, T, np.real(B), np.imag(B))
