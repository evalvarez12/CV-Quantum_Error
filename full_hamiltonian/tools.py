# -*- coding: utf-8 -*-
"""
Utility functions used in CV operations

Created on Tue Feb 19 11:51:17 2019

@author: Eduardo Villasenor 
"""


import numpy as np
import qutip as qt


def tensor(N, operator, pos, Nmodes):
    """
    N:  size of the individual Hilbert space extra dimension
    pos: position of the operator relative to N_list
    
    """
    # The list of extra Hilbert spaces to be tensored
    N_list = [N]*(Nmodes - 1)
    
    
    if N_list[pos:]:
        identity_right = qt.identity(N_list[pos:])
        operator = qt.tensor(operator, identity_right)

    if N_list[:pos]:
        identity_left = qt.identity(N_list[:pos])    
        operator = qt.tensor(identity_left, operator)
    
    return operator
    
def reorder_two_mode_operator(N, op, pos, Nmodes):
        op = qt.tensor([op] + [qt.qeye(N)]*(Nmodes-2))

        # TODO improve this list to permute
        permute_list = list(range(Nmodes))

        permute_list[pos[0]] = 0
        permute_list[0] = pos[0]

        permute_list[pos[1]] = permute_list[1]
        permute_list[1] = pos[1]
        op = op.permute(permute_list)
        return op
    
    
def tmsv_bad_method(N, mean_photon_number):
    r = np.arcsinh(np.sqrt(mean_photon_number))
    psi_theory = qt.tensor(qt.basis(N), qt.basis(N)) * 0

    for i in range(N):
    #    psi_theory += alpha_n(i) * r_nkT(i, 1, T) * qt.tensor(qt.basis(N, i), qt.basis(N, i))
        psi_theory += alpha_n(r, i) * qt.tensor(qt.basis(N, i), qt.basis(N, i))
#        print(i, alpha_n(r, i))
    return psi_theory/psi_theory.norm()
    
    
    
def alpha_n(r, n):
    return (np.tanh(r)**n)/np.cosh(r)