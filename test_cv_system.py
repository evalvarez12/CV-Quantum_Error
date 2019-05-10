# -*- coding: utf-8 -*-
"""
Code to test CV system
Created on Mon Apr 15 10:47:50 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import cv_system as cv
import tools
import matplotlib.pyplot as plt
import unittest



class TestCVSystemsMethods(unittest.TestCase):

   def test_apply_scissor_exact(self):
       N = 3
       sys = cv.System(N, Nmodes=1)
       r = .6
       sys.apply_SMD(r)

       ref_state = sys.state
       k = .05
       sys.apply_scissor_exact(k)

       g = np.sqrt(1/k - 1)
       ref_state = ref_state.data.toarray()
       ref_state[1] = g * ref_state[1]
       ref_state[2:] = 0
       A = np.linalg.norm(ref_state)
       ref_state = ref_state/A
       ref_state = qt.Qobj(ref_state)
       self.assertTrue(sys.state == ref_state)


   def test_apply_scissor_exact2(self):
       N = 3
       sys = cv.System(N, Nmodes=2)
       r = .6
       sys.apply_SMD(r, 1)

       ref_state = sys.state
#       print(ref_state)
       
       k = .05
       sys.apply_scissor_exact(k, 1)
#       print(sys.state)

       g = np.sqrt(1/k - 1)
#       print("Gain:", g)
       ref_state = ref_state.data.toarray()
       ref_state[3] = g * ref_state[3]
       ref_state[6] = 0
       A = np.linalg.norm(ref_state)
       ref_state = ref_state/A
       ref_state = qt.Qobj(ref_state, dims=[[3, 3], [1, 1]])
#       print(ref_state)
       self.assertTrue(sys.state == ref_state)


   def test_apply_loss_channel(self):
       N = 2
       sys = cv.System(N, Nmodes=2)
       r = .6
       sys.apply_SMD(r, 1)

       ref_state = sys.state
       print(ref_state)
       
       eta = 1
       
       sys.apply_loss_channel(eta, 1)
       print(sys.state)
       
       print(ref_state * ref_state.dag())
       
if __name__ == '__main__':
   unittest.main()
