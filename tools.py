# -*- coding: utf-8 -*-
"""
Utility functions used in CV operations

Created on Tue Feb 19 11:51:17 2019

@author: Eduardo Villasenor
"""


import numpy as np
import qutip as qt
from scipy.linalg import block_diag


def tensor(N, operator, pos, Nmodes):
    """
    N:  size of the individual Hilbert space extra dimension
    pos: position of the operator relative to N_list, starts at 0
    Nmodes: number of modes used in the full operator
    """
    if pos >= Nmodes:
        raise ValueError("pos is higher than the number of modes")

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

        permute_list = np.arange(Nmodes)

        # Swap first element
        permute_list[pos[0]] = 0
        permute_list[0] = pos[0]

        # Find where the value 1 is
        ind1 = np.where(permute_list == 1)[0][0]
        # Swap this index with pos[1]
        permute_list[pos[1]] = 1
        permute_list[ind1] = pos[1]
        op = op.permute(permute_list)
        return op


def matrix_sandwich(A, B):
    # Returns A.transpose * B * A
    return np.dot(A.transpose(), np.dot(B, A))


def direct_sum(matrices):
    return block_diag(*matrices)
    

def direct_sum_singles(matrices, positions, Nmodes):
     # Check if arguments make sense
    if max(positions) >= Nmodes:
        raise ValueError("direct_sum: postion out ot numbe blocks")

    identity = np.identity(2)
    positions = np.array(positions)
    all_matrices = [identity] * Nmodes
    for i in range(len(matrices)):
        pos = positions[i]
        all_matrices[pos] = matrices[i]

    return block_diag(*all_matrices)


def reorder_two_mode_symplectic(S, pos, Nmodes):
    # Create first the full S as a list of blocks z
    S_full = np.zeros((2*Nmodes, 2*Nmodes))
    
    # Set diagonal equal to 1
    np.fill_diagonal(S_full, 1)

    p1 = pos[0]*2
    p2 = pos[1]*2

    # Take the blocks of S
    S_full[p1:p1+2, p1:p1+2] = S[0:2, 0:2]
    S_full[p1:p1+2, p2:p2+2] = S[0:2, 2:4]
    S_full[p2:p2+2, p1:p1+2] = S[2:4, 0:2]
    S_full[p2:p2+2, p2:p2+2] = S[2:4, 2:4]
    
    return S_full
    

def kron(matrices):
    if len(matrices) == 2:
        return np.kron(matrices[0], matrices[1])
    else:
        return np.kron(matrices[0], kron(matrices[1:]))


def kron_sum(matrices, positions, Nmodes):
        # Check if arguments make sense
    if max(positions) >= Nmodes:
        raise ValueError("direct_sum: postion out ot numbe blocks")

    N = int(np.sqrt(matrices[0].size))
    identity = np.identity(N)

    full_matrix = 0
    # TODO: remove this for
    for i in range(len(matrices)):
        row = [identity]*Nmodes
        row[positions[i]]  = matrices[i]
        submatrix = kron(row)
        full_matrix += submatrix

    return np.block(full_matrix)
