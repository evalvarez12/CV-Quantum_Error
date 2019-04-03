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