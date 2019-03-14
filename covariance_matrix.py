# -*- coding: utf-8 -*-
"""
Functions to work on CV systems using covariance matrices

Created on Wed Mar 13 11:04:05 2019

@author: Eduardo Villasenor
"""


import numpy as np


#
#a = np.array([[1,2],[2,3]])
#b = a*0
#np.block([[a,b],[b,a]])

def direct_sum(matrices, positions, nblocks):
    # Check if arguments make sense
    if max(positions) >= nblocks:
        raise ValueError("direct_sum: postion out ot numbe blocks")
    
    zero_mat = np.zeros_like(matrices[0])
    N = int(np.sqrt(zero_mat.size))
    identity = np.identity(N)
    
    
    all_rows = []
    # TODO: remove this for
    for i in range(nblocks):
        row = [zero_mat]*nblocks
        if i in positions:
            row[i] = matrices[np.where(i in positions)[0][0]]
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


def beam_splitter_symplectic(theta, positions=[0,1], nmodes=2):
    # The symplectic transformation for  a single mode
    S = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    
    # Apply the direct sum to obtain the full transformation
    S = direct_sum([S, S], positions, nmodes)
    return S

#a = np.array([[1,2],[3,4]])
#
#print(direct_sum([a,a], [0,3], 3))
#
#print(beam_splitter_symplectic(np.pi/4, [0, 1], 2))
    
#w = np.array([[0, 1], [-1, 0]])
#
#print(direct_sum([w, w, w], [0, 1, 2], 3))