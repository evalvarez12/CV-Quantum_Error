# -*- coding: utf-8 -*-
"""
Utility functions used in CV operations

Created on Tue Feb 19 11:51:17 2019

@author: Eduardo Villasenor 
"""


import numpy as np
import qutip as qt


def tensor(operator, N, pos, extraN):
    """
    N:  size of the individual Hilbert space extra dimension
    pos: position of the operator relative to N_list
    
    """
    N_list = [N]*extraN
    
    
    if N_list[:pos]:
        identity_right = qt.identity(N_list[:pos])
        operator = qt.tensor(operator, identity_right)

    if N_list[pos:]:
        identity_left = qt.identity(N_list[pos:])    
        operator = qt.tensor(identity_left, operator)
    
    return operator
    