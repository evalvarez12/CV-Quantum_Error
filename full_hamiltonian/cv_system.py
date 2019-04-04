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
        self.Nmodes = Nmodes


    def add_vacuum(self, Nadd=1):
        self.state = qt.tensor([self.state] + [qt.basis(self.N, 0)]*Nadd)
        self.Nmodes = self.Nmodes + Nadd

    def apply_BS(self, t, pos):
        theta = np.arccos(np.sqrt(t))
        U = bs.beam_splitter_U(self.N, theta, pos, self.Nmodes)

        if self.state.isket:
            self.state = U * self.state
        else:
            self.state = U * self.state * U.dag()
            
    
    def save_state(self):
        self.state_saved = self.state
        
    
    def load_state(self):
        if self.state_saved is None:
            raise AttributeError("There is not a saved state")
        self.state = self.state_saved
        
        
    def load_state_del(self):
        if self.state_saved is None:
            raise AttributeError("There is not a saved state")
        self.state = self.state_saved
        del self.state_saved
        
    def apply_TMS(self, mphoton, pos):
        r = np.arcsinh(np.sqrt(mphoton))
        S = ops.tmsqueeze(self.N, r, pos, self.Nmodes)
        
        if self.state.isket:
            self.state = S * self.state
        else:
            self.state = S * self.state * S.dag()
            
    
    def apply_SMS(self, mphoton, pos):
        r = np.arcsinh(np.sqrt(mphoton))
        S = ops.tmsqueeze(self.N, r, pos, self.Nmodes)
        
        if self.state.isket:
            self.state = S * self.state
        else:
            self.state = S * self.state * S.dag()
    
    
    def add_TMSV(self, mphoton):
        r = np.arcsinh(np.sqrt(mphoton)) 
        state_aux = qt.tensor(qt.basis(self.N), qt.basis(self.N))
        S = ops.tmsqueeze(self.N, r)
        
        state_aux = S * state_aux 
        if not self.state.isket:
            state_aux = state_aux * state_aux.dag()
        self.state = qt.tensor(self.state, state_aux)
        self.Nmodes = self.Nmodes + 2

    
    def collapse_fock_state(self, fock, pos):
        P = qt.basis(self.N, fock).dag()
        P = tools.tensor(self.N, P, pos, self.Nmodes)

        # Compute the probability of this collapse when measured
        if self.state.isket:
            p = (P * self.state).norm()
        else:
            # TODO: Check this
            p = (P.dag() * P * self.state).tr()

        self.state = (P * self.state)/p
        self.Nmodes = self.Nmodes - 1
        return p


    def get_simple_CM_V(self, mode):
        a = qt.destroy(self.N)
        a = tools.tensor(self.N, a, mode, self.Nmodes)

        v = 1 + 2 * self.state.dag() * a.dag() * a * self.state
        return v


    def get_simple_CM_C(self, modes):
        a = qt.destroy(self.N)
        a1 = tools.tensor(self.N, a, modes[0], self.Nmodes)
        a2 = tools.tensor(self.N, a, modes[1], self.Nmodes)

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

        # TODO: check this factor of 2
        cm = 2 * qt.covariance_matrix(basis, self.state)
        return cm
