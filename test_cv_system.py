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
       N = 10
       sys = cv.System(N, Nmodes=2)
       r = .6
       sys.apply_SMD(r, 1)

       ref_state = sys.state
#       print(ref_state)
       
       k = .05
       sys.apply_scissor_exact(k, 1)
#       print(sys.state)

       # Set up reference state
       g = np.sqrt(1/k - 1)
#       print("Gain:", g)
       ref_state = ref_state.data.toarray()
       ref_state[N] = g * ref_state[N]
       
       mask = np.ones_like(ref_state, dtype=bool)
       mask[0] = False
       mask[N] = False
       ref_state[mask] = 0
       A = np.linalg.norm(ref_state)
       ref_state = ref_state/A
       ref_state = qt.Qobj(ref_state, dims=[[N, N], [1, 1]])
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
       
       
   def test_tmsv_to_smsvs(self):
       N = 40
       mpn = 1
       t = np.pi/4
        
       r = np.arcsinh(np.sqrt(mpn))
        
       ## Initialize TMSV state
       sys = cv.System(N, Nmodes=2, cm=True)
       sys.apply_TMS(r, [0, 1])
                      
       sys.set_quadratures_basis()
       CM = sys.get_full_CM()
       np.testing.assert_almost_equal(sys.cm, CM, decimal=4)
        
       # Apply BS operation
       sys.apply_BS(t, [0, 1])
        
#       sys.set_quadratures_basis()
       CM = sys.get_full_CM()
       np.testing.assert_almost_equal(sys.cm, CM, decimal=4)
        
       # Second BS to revert to TMSV
       sys.apply_BS(t, [0, 1])
        
       sys.set_quadratures_basis()
       CM = sys.get_full_CM()
       np.testing.assert_almost_equal(sys.cm, CM, decimal=4)
       
       
   def test_smsvs_to_tmsv(self):
       N = 40
       mpn = 1
       t = np.pi/4
        
       r = np.arcsinh(np.sqrt(mpn))       
       sys = cv.System(N, Nmodes=2, cm=True)
       #z = mpn*np.exo(1j*np.pi/4)
       sys.apply_SMS(-r, 0)
       sys.apply_SMS(r, 1)
        
       sys.set_quadratures_basis()
       CM = sys.get_full_CM()
       np.testing.assert_almost_equal(sys.cm, CM, decimal=4)
        
       sys.apply_BS(t, [0, 1])
       CM = sys.get_full_CM()
       np.testing.assert_almost_equal(sys.cm, CM, decimal=4)
        
       sys.apply_BS(t, [0, 1])
       sys.set_quadratures_basis()
       CM = sys.get_full_CM()
       np.testing.assert_almost_equal(sys.cm, CM, decimal=4)
if __name__ == '__main__':
   unittest.main()
