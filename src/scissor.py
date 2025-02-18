# -*- coding: utf-8 -*-
"""
Operations to realize entanglement distillation via quantum scisors

Created on Thu Feb 28 11:04:31 2019

@author: Eduardo Villasenor
"""

import numpy as np
import qutip as qt
from . import operations as ops
from . import beam_splitter as bs
from . import tools







def exact(rho_in, kappa):
    """
             D
             |    
  c=in  ->---/-----D
             |
             |
  b=|0> ->---/-----> out
             |
             ^
            a=|1>
            
    H = in * in_extra : only applies to the leftmost mode of the in state 
    """
    N = rho_in.dims[0][0]
    N_modes = len(rho_in.dims[0])
        
    theta1 = np.arccos(np.sqrt(kappa))
    theta2 = np.pi/4
    
    state_a = qt.basis(N, 1) * qt.basis(N, 1).dag()
    state_b = qt.basis(N) * qt.basis(N).dag()
    
    rho = qt.tensor([state_a, state_b, rho_in])
#    print(rho, rho.tr())
    rho = bs.tritter_applyU(rho, theta1, theta2)
#    print(rho, rho.tr())
    
    # Define the proyectors, in this case to |10>    
    projector0 = qt.basis(N, 0).dag()
    projector1 = qt.basis(N, 1).dag()
    projector = qt.tensor([projector0, qt.identity(N), projector1])
    if N_modes > 1:
        projector = tools.tensor(projector, N, 0, N_modes)

#    print(projector)    
    rho = projector * rho * projector.dag()
    rho = rho/rho.tr()
    
    return rho

def NLA(rho_in, kappa, mu_aux):
    """
             D
             |    
  c=in  ->---/----D
             |
             |
  b=|0> ->---/----> out
             |
             ^
           a=TMSV ->----D
            
    H = in * in_extra : only applies to the leftmost mode of the in state 
    """
    N = rho_in.dims[0][0]
    N_modes = len(rho_in.dims[0])
    
    theta1 = np.arccos(np.sqrt(kappa))
    theta2 = np.pi/4
    
    state_b = qt.basis(N) * qt.basis(N).dag()
    
    
    state_a = qt.tensor(state_b, state_b)
    S = ops.tmsqueeze(N, mu_aux)
    state_a = S * state_a * S.dag()
#    state_a = qt.basis(N, 1) * qt.basis(N, 1).dag()
    
    rho = qt.tensor([state_a, state_b, rho_in])
#    print(rho, rho.tr())
    rho = bs.tritter_applyU(rho, theta1, theta2, pos=[0, 2, 3])
#    print(rho, rho.tr())
    
    # Define the proyectors, in this case to |10>    
    projector_0 = qt.basis(N, 0).dag()
    projector_on = ops.photon_on_projector(N)
#    proyector10 = qt.tensor([proyector0, qt.identity(N), proyector1])
    projector = qt.tensor([projector_0, projector_on, qt.identity(N), projector_on])
    
    if N_modes > 1:
        projector = tools.tensor(projector, N, 0, N_modes)
    
    rho = projector * rho * projector.dag()
    rho = rho/rho.tr()
    
    return rho


