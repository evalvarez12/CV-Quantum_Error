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
       N = 10
       sys = cv.System(N, Nmodes=1)
       r = .6
       sys.apply_SMD(r)

       ref_state = sys.state
       # print(ref_state)
       k = .05
       sys.apply_scissor_exact(k)
       # print(sys.state)

       g = np.sqrt(1/k - 1)
       # print("Gain:", g)
       ref_state = ref_state.data.toarray()
       ref_state[1] = g * ref_state[1]
       ref_state[2:] = 0
       A = np.linalg.norm(ref_state)
       ref_state = ref_state/A
       # print(ref_state)
       ref_state = qt.Qobj(ref_state)
       self.assertTrue(sys.state == ref_state)


if __name__ == '__main__':
   unittest.main()
