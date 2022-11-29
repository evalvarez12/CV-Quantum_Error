#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 12:06:37 2022

@author: eduardo
"""

import numpy as np

I = np.identity(2)
Z = np.diag([1, -1])
O = np.zeros((2,2))

def BS_mat(T, pos1=0, pos2=1, n=2):
    pos1 *= 2
    pos2 *= 2
    bs = np.identity(2*n)
    bs[pos1:pos1+2, pos1:pos1+2] = np.sqrt(T)*I
    bs[pos1:pos1+2, pos2:pos2+2] = np.sqrt(1-T)*I
    bs[pos2:pos2+2, pos1:pos1+2] = -np.sqrt(1-T)*I
    bs[pos2:pos2+2, pos2:pos2+2] = np.sqrt(T)*I
    return bs

def swap(A, pos1, pos2):
    As = np.copy(A)
    pos1 *= 2
    pos2 *= 2
    As[pos1:pos1+2, pos1:pos1+2] = A[pos2:pos2+2, pos2:pos2+2]
    As[pos1:pos1+2, pos2:pos2+2] = A[pos2:pos2+2, pos1:pos1+2]
    As[pos2:pos2+2, pos1:pos1+2] = A[pos1:pos1+2, pos2:pos2+2]
    As[pos2:pos2+2, pos2:pos2+2] = A[pos1:pos1+2, pos1:pos1+2]
    return As


def TMSV_mat(v):
    return np.block([[v*I, np.sqrt(v**2-1)*Z],
                     [np.sqrt(v**2-1)*Z, v*I]])

def dir_sum2(A, B):
    dsum = np.zeros(np.add(A.shape,B.shape))
    dsum[:A.shape[0],:A.shape[1]]=A
    dsum[A.shape[0]:,A.shape[1]:]=B
    return dsum


def dir_sum(As):
   if len(As) == 1:
       return As[0]
   else:
       ds = dir_sum2(As[0],As[1])
       return dir_sum([ds] + As[2:])
        
"""
SYSTEM
               
A -----/-------/------/------ B
       J      E1      F
       K      E2      G
"""


ta = 0.8
tb = 0.8
te = 0.8
va = 1.2
vb = 1
ve = 1
v = 1.4

Sb = BS_mat(tb)
Sa = BS_mat(ta)
Se = BS_mat(te)
JK = TMSV_mat(va)
FG = TMSV_mat(vb)
EE = TMSV_mat(ve)
AB = TMSV_mat(v)


S1 = dir_sum([AB, JK])

BS1 = dir_sum([I,  Sa, I]) 
S2 = np.dot(np.dot(BS1.transpose(), S1), BS1)

S3 = dir_sum([S2, EE])
BS2 = BS_mat(te, 1, 4, 6)
S4 = np.dot(np.dot(BS2.transpose(), S3), BS2)



# S2 = S2[2]

# # Eve's modes
# S3 = dir_sum([S2, EE])
# BS2 = dir_sum([I,I,Se,I])
