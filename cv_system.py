# -*- coding: utf-8 -*-
"""
Class to simulate a CV quantum system

Created on Wed Apr  3 09:29:55 2019

@author: Eduardo Villasenor
"""

import operations as ops
import symplectic as sym
import beam_splitter as bs
import tools
import theory
import qutip as qt
import numpy as np



class System:
    def __init__(self, N, Nmodes=0, cm=False, cm_only=False):
        self.N = N
        if Nmodes != 0 and not cm_only:
            self.state = qt.tensor([qt.basis(self.N, 0)]*Nmodes)
            self.Nmodes = Nmodes

        if cm:
            eye = np.eye(2)
            self.cm = tools.direct_sum([eye]*Nmodes)
        else:
            self.cm = None

        self.quad_basis = None

    def add_vacuum(self, Nadd=1):
        state_add = qt.tensor([qt.basis(self.N, 0)]*Nadd)
        if not self.state.isket:
            state_add = state_add * state_add.dag()
            
        self.state = qt.tensor(state_add, self.state)
        self.Nmodes = self.Nmodes + Nadd

    
    def add_state(self, state):
        N_add = len(state.dims[0])
        self.state = qt.tensor(state, self.state)
        self.Nmodes += N_add

    def apply_BS(self, z, pos=[0,1]):
#        theta = np.arccos(np.sqrt(t))
#        U = bs.beam_splitter_U(self.N, theta, pos, self.Nmodes)
        U = ops.beam_splitter(self.N, z, pos, self.Nmodes)

        if self.state.isket:
            self.state = U * self.state
        else:
            self.state = U * self.state * U.dag()

        if self.cm is not None:
            S = sym.beam_splitter(z, pos, self.Nmodes)
            self.cm = tools.matrix_sandwich(S, self.cm)


    def apply_TMS(self, r, pos=[0,1]):
        U = ops.tmsqueeze(self.N, r, pos, self.Nmodes)

        if self.state.isket:
            self.state = U * self.state
        else:
            self.state = U * self.state * U.dag()

        if self.cm is not None:
            S = sym.two_mode_squeeze(r, pos, self.Nmodes)
            self.cm = tools.matrix_sandwich(S, self.cm)


    def apply_SMD(self, alpha, pos=0):
        U = qt.displace(self.N, alpha)
        U = tools.tensor(self.N, U, pos, self.Nmodes)

        if self.state.isket:
            self.state = U * self.state
        else:
            self.state = U * self.state * U.dag()


    def apply_SMS(self, r, pos=0):
        U = ops.squeeze(self.N, r, pos, self.Nmodes)

        if self.state.isket:
            self.state = U * self.state
        else:
            self.state = U * self.state * U.dag()

        if self.cm is not None:
            S = sym.single_mode_squeeze(r, pos, self.Nmodes)
            self.cm = tools.matrix_sandwich(S, self.cm)


    def apply_loss_channel(self, eta, pos):
        eta_theta = np.arccos(np.sqrt(eta))
        self.add_vacuum()
        self.apply_BS(eta_theta, [pos, self.Nmodes-1])
        
        # Trace out loss channel mode - NOTE in qutip the last one is index 0
        self.state = self.state.ptrace(range(1, self.Nmodes))
        self.Nmodes = self.Nmodes - 1

    def add_TMSV(self, r):
        state_aux = qt.tensor(qt.basis(self.N), qt.basis(self.N))
        S = ops.tmsqueeze(self.N, r)

        state_aux = S * state_aux
        if not self.state.isket:
            state_aux = state_aux * state_aux.dag()
        self.state = qt.tensor(state_aux, self.state)
        self.Nmodes = self.Nmodes + 2

        if self.cm is not None:
            S = sym.two_mode_squeeze(r)
            cm_add = np.dot(S.transpose(), S)
            self.cm = tools.direct_sum([self.cm, cm_add])




    def collapse_fock_state(self, fock, pos):
        P = qt.basis(self.N, fock).dag()
        P = tools.tensor(self.N, P, pos, self.Nmodes)

        # Compute the probability of this collapse when measured
        if self.state.isket:
            p = (P * self.state).norm()
        else:
            # TODO: Check this
            p = (P.dag() * P * self.state).tr()

#        if p == 0:
#            p = 1
        self.state = (P * self.state)/p
        self.Nmodes = self.Nmodes - 1
        return p


    def apply_scissor_exact(self, t, pos=0):
        # Tritter parameters
        theta1 = np.pi/4
        theta2 = np.arccos(np.sqrt(t))

        # Add extra state |10>
        extra_psi = qt.tensor(qt.basis(self.N, 1), qt.basis(self.N, 0))
        if not self.state.isket:
            extra_psi = extra_psi * extra_psi.dag()
        self.state = qt.tensor(extra_psi, self.state)
        Nmodes = self.Nmodes + 2

        # Apply tritter operator
        tritter_pos=[Nmodes-1, Nmodes-2, pos]
        
        ######
#        tritter_pos.reverse()
        
        U = ops.tritter(self.N, theta1, theta2, tritter_pos, Nmodes)
        # print(Nmodes, theta1, theta2, U)
        if self.state.isket:
            self.state = U * self.state
        else:
            self.state = U * self.state * U.dag()

        # Define the proyectors, in this case to |10>
        projector0 = qt.basis(self.N, 0).dag()
        projector1 = qt.basis(self.N, 1).dag()
        collapse_pos = [pos, Nmodes-1]
        projector = tools.tensor_singles(self.N, [projector0, projector1], collapse_pos, Nmodes)

        if self.state.isket:
            self.state = projector * self.state
        else:
            self.state = projector * self.state * projector.dag()

        # Normalize the state
        if self.state.isket:
            p_success = self.state.norm()
        else:
            p_success = self.state.tr()

        if p_success == 0:
            p_success = 1
        self.state = self.state/p_success
        
        # Permute the state to return it to its original ordering
        p_list = np.arange(self.Nmodes)
        p_list[pos] = self.Nmodes - 1
        p_list[-1] = pos
        self.state = self.state.permute(p_list)
        return p_success
    
    
    def apply_scissor(self, t, r_aux, pos=0):
        # Tritter parameters
        theta1 = np.pi/4
        theta2 = np.arccos(np.sqrt(t))

        # Add extra vacuum and TMSV states
        self.add_TMSV(r_aux)
        self.add_vacuum()
        
        # Apply tritter operator
        tritter_pos=[self.Nmodes-1, self.Nmodes-2, pos]
        U = ops.tritter(self.N, theta1, theta2, tritter_pos, self.Nmodes)
        # print(Nmodes, theta1, theta2, U)
        if self.state.isket:
            self.state = U * self.state
        else:
            self.state = U * self.state * U.dag()

        # Define the proyectors, in this case to |10>
        projectorOFF = qt.basis(self.N, 0).dag()
        projectorON = ops.photon_on_projector(self.N)
        collapse_pos = [pos, self.Nmodes-1, self.Nmodes-3]
        projector = tools.tensor_singles(self.N, [projectorOFF, projectorON, projectorON], collapse_pos, self.Nmodes)

        if self.state.isket:
            self.state = projector * self.state
        else:
            self.state = projector * self.state * projector.dag()

        # Return Nmodes to original value
        self.Nmodes = self.Nmodes - 3
        
        # Normalize the state
        if self.state.isket:
            p_success = self.state.norm()
        else:
            p_success = self.state.tr()

        if p_success == 0:
            p_success = 1
        self.state = self.state/p_success
        
        # Permute the state to return it to its original ordering
        p_list = np.arange(self.Nmodes)
        p_list[pos] = self.Nmodes - 1
        p_list[-1] = pos
        self.state = self.state.permute(p_list)
        return p_success


    def ptrace(self, pos_keep):
        self.Nmodes = len(pos_keep)
        self.state = self.state.ptrace(pos_keep)
        

    def get_simple_CM_V(self, mode):
        a = qt.destroy(self.N)
        a = tools.tensor(self.N, a, mode, self.Nmodes)

        v = 1 + 2 * self.state.dag() * a.dag() * a * self.state
        return v.norm()


    def get_simple_CM_C(self, modes):
        a = qt.destroy(self.N)
        a1 = tools.tensor(self.N, a, modes[0], self.Nmodes)
        a2 = tools.tensor(self.N, a, modes[1], self.Nmodes)

        c = self.state.dag() * (a1*a2 + a1.dag()*a2.dag()) * self.state
        return c.norm()


    def set_quadratures_basis(self):
        # basis = {x1, p1, x2, p2, ....., xn, pn}
        a = qt.destroy(self.N)
        x = (a + a.dag())/np.sqrt(2)
        p = 1j*(a - a.dag())/np.sqrt(2)

        basis = []
        # TODO: check the order of the basis
        for i in range(self.Nmodes):
            basis += [tools.tensor(self.N, x, i, self.Nmodes),
                      tools.tensor(self.N, p, i, self.Nmodes)]

        self.quad_basis = basis


    def get_CM_entry(self, indeces):
        am = self.quad_basis[indeces[0]]
        an = self.quad_basis[indeces[1]]
        a = qt.expect(am * an +  an * am, self.state)/2
        b = qt.expect(am, self.state) * qt.expect(an, self.state)
        Vmn = a - b
        return 2 * Vmn


    def get_full_CM(self):
        if self.quad_basis is None:
            self.set_quadratures_basis()
        
        # TODO: check this factor of 2
        cm = 2 * qt.covariance_matrix(self.quad_basis, self.state)
        return cm.astype(float)


    def check_CM(self):
        cm = self.get_full_CM()
        return np.array_equal(cm, self.cm)


    def replace_current_state_w_bad_TMSV(self, mean_photon_number):
        self.state = tools.tmsv_bad_method(self.N, mean_photon_number)
        self.Nmodes = 2


    def add_bad_TMSV(self, mean_photon_number):
        # NOTE: just for comparisson - shuold delete
        self.state = qt.tensor(self.state, theory.tmsv_state(self.N, mean_photon_number))
        self.Nmodes += 2


    def save_state(self):
        self.saved = [self.state, self.Nmodes, self.cm]


    def load_state(self):
        if self.saved is None:
            raise AttributeError("There is not a saved state")
        self.state = self.saved[0]
        self.Nmodes = self.saved[1]

        if self.cm is not None:
            self.cm = self.saved[2]


    def load_state_del(self):
        self.load_state()
        del self.saved
        
        
    def reset_state(self):
        self.state = qt.tensor([qt.basis(self.N, 0)]*self.Nmodes)

        if self.cm is not None:
            eye = np.eye(2)
            self.cm = tools.direct_sum([eye]*self.Nmodes)
