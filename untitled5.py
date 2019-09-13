# -*- coding: utf-8 -*-
"""
Random code to test the purification of covariance matrix
Created on Tue Sep 10 17:41:25 2019

@author: z5239621
"""

import numpy as np

alpha = np.sqrt(5)

a = (1 + 2*alpha)
b = 2*(np.sqrt(alpha**2 + alpha))

Z = np.array([[1, 0], [0, -1]])
I = np.identity(2)

cm_TMSV = np.block([[a*I, b*Z],[b*Z, a*I]])
print("CM_TMSV:", np.linalg.det(cm_TMSV))

ts = 0.9
t = alpha/(1+alpha)*ts
N = 1

x = 1 + 2*(N+t)/(1-t)
y = 1 + 2*((N+1)*t)/(1-t)
z = 2*np.sqrt(t)*(N+1)/(1-t)

cm_PSS = np.block([[x*I, z*Z],[z*Z, y*I]])
print("CM_PSS:", np.linalg.det(cm_PSS))

te = 1
v = 1.1
e = 0.1
nr = 0.95
nd = 0.68

z2 = z * np.sqrt(te)
y2 = te*(y + e) + (1-te)

cag = -np.sqrt(1-nd)*z2
cgh = np.sqrt(nd*(v**2-1))

gb3 = (nd*y2 + (1-nd)*v)*I

gagh = np.block([[x*I, cag*Z, 0*I],[cag*Z, (nd*v + (1-nd)*y2)*I, cgh*Z], [0*I, cgh*Z, v*I]])
saghb3 = np.block([np.sqrt(nd)*z2*Z, np.sqrt((1-nd)*nd)*(v-y2)*Z, np.sqrt((1-nd)*(v**2 - 1)) * Z ]).transpose()

yaghb3 = np.block([[gagh, saghb3],[saghb3.transpose(), gb3]])
print("CM_purified:", np.linalg.det(yaghb3))