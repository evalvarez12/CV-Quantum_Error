# -*- coding: utf-8 -*-
"""
Tools to assist on covariance matrix operations

Created on Fri Mar 15 12:14:36 2019

@author: Eduardo Villasenor
"""

import numpy as np

def direct_sum(matrices, positions, nblocks):
    # Check if arguments make sense
    if max(positions) >= nblocks:
        raise ValueError("direct_sum: postion out ot numbe blocks")
    
    zero_mat = np.zeros_like(matrices[0])
    N = int(np.sqrt(zero_mat.size))
    identity = np.identity(N)
    
    positions = np.array(positions)
    
    all_rows = []
    # TODO: remove this for
    for i in range(nblocks):
        row = [zero_mat]*nblocks
        if i in positions:
            pos = np.where(i == positions)[0][0]
            row[i] = matrices[pos]
        else:
            row[i] = identity
        
        all_rows += [row]
        
    return np.block(all_rows)
        

def kron(matrices):
    if len(matrices) == 2:
        return np.kron(matrices[0], matrices[1])
    else:
        return np.kron(matrices[0], kron(matrices[1:]))

def kron_sum(matrices, positions, nblocks):
        # Check if arguments make sense
    if max(positions) >= nblocks:
        raise ValueError("direct_sum: postion out ot numbe blocks")
    
    N = int(np.sqrt(matrices[0].size))
    identity = np.identity(N)
    
    full_matrix = 0
    # TODO: remove this for
    for i in range(len(matrices)):
        row = [identity]*nblocks
        row[positions[i]]  = matrices[i]
        submatrix = kron(row)
        full_matrix += submatrix
        
    return np.block(full_matrix)
