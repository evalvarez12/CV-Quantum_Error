# -*- coding: utf-8 -*-
"""
Test routines for quantum scissors

Created on Fri Mar  8 10:58:55 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import scissor
sys.path.append("..") 



N =5
kappa = .05

a = qt.displace(N,1)

vacuum = qt.basis(N) * qt.basis(N).dag()
D = qt.displace(N, 1)

rho = D * vacuum * D.dag()

#state = qt.rand_dm(N)

#rho = qt.tensor(rho, vacuum)
print(rho) 

result = scissor.exact(rho, kappa)

print(result)

mu_aux = 0.1
rho1NLA = scissor.NLA(rho, kappa, mu_aux)
print(rho1NLA)

