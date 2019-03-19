# -*- coding: utf-8 -*-
"""
Recreating the results of GRCI in quantum scissor paper

Created on Mon Mar 18 15:24:59 2019

@author: Eduardo Villasenor
"""


import numpy as np
import covariance_matrix as cm
import tools_cm as ts



# Define parameters
mu = 1.
kappa = 0.5

eta = 0.01
mu_aux = 0.01
identity = np.identity(2)


cm_ab = cm.tmsv_cm(mu)
cm_c = identity
cm_de = cm.tmsv_cm(mu_aux)


cm_0 = np.kron(np.kron(cm_ab, cm_c), cm_de)

cm_1 = cm.pure_loss(cm_0, 1, 5, eta) 

bs_bd = cm.beam_splitter_symplectic(0.5, [1, 3], 5)
bs_cd = cm.beam_splitter_symplectic(kappa, [2, 3], 5)

bs = np.dot(bs_bd, bs_cd)
bst = np.dot(np.transpose(bs_cd), np.transpose(bs_bd))

cm_2 = np.dot(np.dot(bs, cm_1), bst)