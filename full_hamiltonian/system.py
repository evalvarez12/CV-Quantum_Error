# -*- coding: utf-8 -*-
"""
Class to simulate a CV quantum system

Created on Wed Apr  3 09:29:55 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import operations as ops
import beam_splitter as bs
import tools


class System:
    def __init__(self, N, Nmodes):
        self.N = N
        self.state = qt.tensor([qt.basis(self.N, 0)]*Nmodes)
        self.ket = True
        self.Nmodes = Nmodes


    def add_vacuum(self, Nadd):
        self.state = qt.tensor([self.state] + [qt.basis(self.N, 0)*Nadd])
        self.Nmodes = self.Nmodes + Nadd

    def apply_BS(self, theta, pos):
        U = bs.beam_splitter_U(self.N, theta, pos, self.Nmodes)

        if self.ket:
            self.state = U * self.state
        else:
            self.state = U * self.state * U.dag()

    def collapse_fock_state(self, fock, pos):
        P = qt.basis(self.N, fock).dag()
        P = tools.tensor(self.N, P, pos, self.Nmodes)

        # Compute the probability of this collapse when measured
        if self.ket:
            p = (P * self.state).norm()
        else:
            # TODO: Check this
            p = (P.dag() * P self.state).tr()

        self.state = (P * self.state)/p
        return p


    def get_simple_CM_V(self, mode):
        a = qt.destroy(self.N)
        a = tools.tensor(self.N, a, mode, self.Nmodes)

        v = 1 + 2 * self.state.dag() * a.dag() * a * self.state
        return v


    def get_simple_CM_C(self, modes):
        a = qt.destroy(self.N)
        a1 = tools.tensor(self.N, a, modes[0], self.Nmodes)
        a2 = tools.tensor(self.N, a, modes[0], self.Nmodes)

        c = self.state.dag() * (a1*a2 + a1.dag()*a2.dag()) * self.state
        return c


    def get_full_CM(self):
        a = qt.destroy(self.N)
        x = (a + a.dag())/np.sqrt(2)
        p = 1j*(-a + a.dag())/np.sqrt(2)

        basis = []
        for i in range(self.Nmodes):
            basis += [tools.tensor(self.N, x, i, self.Nmodes),
                      tools.tensor(self.N, p, i, self.Nmodes)]

        cm = qt.covariance_matrix(basis, self.state)
        return cm
