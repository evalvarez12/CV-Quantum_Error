# -*- coding: utf-8 -*-
"""
Routines to play with covariance matrices
Created on Tue Apr  2 13:35:01 2019

@author: Eduardo Villasenor 
"""

import qutip as qt
import numpy as np

N = 20

vacuum = qt.basis(N, 0)

a = qt.create(N).dag()

alpha = 2 + 1j
D = qt.displace(N, alpha)

x = (a + a.dag())/np.sqrt(2)
p = 1j*(-a + a.dag())/np.sqrt(2)

basis = [x,p] 

state = D * vacuum

print("<a>:", qt.expect(a, state))
print("<x>:", qt.expect(x, state))
print("<p>:", qt.expect(p, state))


cm = qt.covariance_matrix(basis, state)
print("CM:", cm)

print("-----------------------------------------------------")

N = 50

vacuum = qt.basis(N, 0)
vacuum = qt.tensor(vacuum, vacuum)

a1 = qt.tensor(qt.create(N), qt.qeye(N)).dag()
a2 = qt.tensor(qt.qeye(N), qt.create(N)).dag() 

eps = 2
S = qt.squeezing(a1, a2, eps)

x1 = (a1 + a1.dag())/np.sqrt(2)
p1 = 1j*(-a1 + a1.dag())/np.sqrt(2)

x2 = (a2 + a2.dag())/np.sqrt(2)
p2 = 1j*(-a2 + a2.dag())/np.sqrt(2)


basis = [x1, p1, x2, p2] 

state = S * vacuum
#state = qt.tensor(qt.rand_dm(N), qt.rand_dm(N))

#print("<a>:", qt.expect(a, state))
#print("<x>:", qt.expect(x, state))
#print("<p>:", qt.expect(p, state))


cm = qt.covariance_matrix(basis, state)
print("CM:", cm)


print((qt.expect(2*a1.dag()*a1, state) + 1)/2)
print(qt.expect(a1*a2 + a1.dag()*a2.dag(), state)/2)
