# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 16:13:07 2019

@author: z5239621
"""
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import scipy.linalg as la

######################################### EB
def abc_EB(V, t, e):
    # Vmod to V
    V = V + 1

    a = V
    b = t*(V-1) + 1 + e
    c = np.sqrt(t*(V**2 - 1))
    return a, b, c


def I_shared_EB(V, t, e, m=1):
    a, b, c = abc_EB(V, t, e)

    if m == 1:
        Vba = b - c**2/a
    elif m == 2:
        b = (b + m - 1)/m
        Vba = b - c**2/(a)

    I = .5 * np.log2(b / Vba)
    return I


def Holevo_EB(V, t, e , m=1):
    a, b, c = abc_EB(V, t, e)


    z = np.sqrt((a+b)**2 - 4*c**2)
    v1 = .5 * (z + (b - a))
    v2 = .5 * (z - (b - a))

    if m == 1:
        v3 = np.sqrt(a * (a - c**2/b))
    elif  m == 2:
        v3 = a - c**2 / (b+1)

#    print(v1, v2, v3)
#    eign_check(a,b,c,m)
    return g(v1) + g(v2) - g(v3)


def key_rate_EB(V, t, e, m=1):
    b = 0.95
    I = I_shared_EB(V, t, e, m)
    X = Holevo_EB(V, t, e, m)
    kr = b*I - X

    print("EB I:", I)
    print("EB X:", X)

    return kr


################################################ PM

def abc_PM(V, t, e, m):
    a, b, c = abc_EB(V, t, e)

    r = np.sqrt(m *(V + 2)/V)
    c = c/r
    a = a - 1
    return a, b, c


def I_shared_PM(V, t, e, m=1):
    I = (m/2) * np.log2(1 + (((t/m) * V ) / (1 + e/m)))
#    I =  np.log2(1 + ((t * V )/ 2* (1 + e/2)))

    return I


def I_shared_PM2(V, t, e, m=1):
    a, b, c = abc_PM(V, t, e, m)

    if m == 1:
        Vba = b - c**2/a
    elif m == 2:
        b = (b + m - 1)/m
        Vba = b - c**2/(a)

    I = m/2 * np.log2(b / Vba)
    return I


def Holevo_PM(V, t, e , m=1):
    a, b, c = abc_EB(V, t, e)

    z = np.sqrt((a+b)**2 - 4*c**2)
    v1 = .5 * (z + (b - a))
    v2 = .5 * (z - (b - a))
    if m == 1:
        v3 = np.sqrt(a * (a - c**2/b))
    elif  m == 2:
        v3 = a - c**2 / (b+1)

#    print(v1, v2, v3)
#    eign_check(a,b,c,m)

    return g(v1) + g(v2) - g(v3)


def key_rate_PM(V, t, e, m=1):
    b = 0.95
    I = I_shared_PM(V, t, e, m)

    # I2 = I_shared_PM2(V, t, e, m)
    # print("Compare I", I, I2)

    X = Holevo_PM(V, t, e, m)

#    print("PM I:", I)
#    print("PM X:", X)

    kr = b*I - X
    return kr

def key_rate_PM_beta(V, t, e, b, m=1):
    I = I_shared_PM(V, t, e, m)
    X = Holevo_PM(V, t, e, m)
    kr = b*I - X
    return kr


####################################################

def eign_check(a, b, c, m=1):
    I = np.eye(2)
    S = np.diag([1, -1])

    Mab = np.block([[I*a, S*c],[S*c, I*b]])
    if m == 1:
        Ma_b = a*I - np.diag([1/b, 0]) * c**2
    elif m == 2:
        Ma_b =  a*I - np.diag([1, 1]) * (c**2/(b+1))


    v = symplectic_eigenvalues(Mab)
    v2 = symplectic_eigenvalues(Ma_b)
    print(v)
    print(v2)


def symplectic_eigenvalues(cm):
    N = int(cm.shape[0]/2)
    omega = symplectic_form(N)
    eigvals = np.linalg.eigvals(1j*np.dot(omega, cm)).real
    eigvals = eigvals[eigvals > 0]
    return eigvals


def symplectic_form(N):
    w = np.array([[0, 1], [-1, 0]])

    # Use direct sum to calculate the symplectic form
    return la.block_diag(*[w]*N)


def g(v):
    return ((v+1)/2)*(np.log2((v+1)/2)) - ((v-1)/2)*(safe_log2((v-1)/2))

def safe_log2(x, minval=0.0000000001):
    return np.log2(x.clip(min=minval))


def avg_opt_key_rate(t, e, m=1):
    t_avg = np.average(t)
    V_opt, _ = opt_key_rate(np.array([t_avg]), e, m)
    return V_opt, key_rate_PM(V_opt, t, e, m)


def opt_key_rate(t, e, m=1):
    V_max = []
    for i in t:
        K = lambda V : -key_rate_PM(V, i, e, m)
        res = op.minimize(K, .5)
        V_max += [res['x'][0]]
#        print(res)
    V_max = np.array(V_max)
    K = key_rate_PM(V_max, t, e, m)
    return V_max, K


def opt_key_rate_eff(t, tsqrt, e, m=1):
    Vart = t - tsqrt**2
    K = lambda V : -key_rate_PM(V, tsqrt**2, Vart*V + e, m)
    res = op.minimize(K, .5)
    V_opt = res['x'][0]
#        print(res)
    K = key_rate_PM(V_opt, tsqrt**2, Vart*V_opt + e, m)
    return V_opt, K


def opt_key_rate_e(t, e, m=1):
    V_max = []
    for i in e:
        K = lambda V : -key_rate_PM(V, t, i, m)
        res = op.minimize(K, .5)
        V_max += [res['x'][0]]
#        print(res)
    V_max = np.array(V_max)
    K = key_rate_PM(V_max, t, e, m)
    return V_max, K

def opt_key_rate_b(t, e, b, m=1):
    V_max = []
    for i in b:
        K = lambda V : -key_rate_PM_beta(V, t, e, i, m)
        res = op.minimize(K, .5)
        V_max += [res['x'][0]]
#        print(res)
    V_max = np.array(V_max)
    K = key_rate_PM_beta(V_max, t, e, b, m)
    return V_max, K


#####################################################################

#
#
##################### CHECK EB
#t = 0.1
#e = 0.01
#V = 1
#print("Key rate EB:", key_rate_EB(V, t, e))
#
#
##x = np.linspace(0.001, 50, 5000)
##t = .1
##for e in [0.0, 0.01, 0.02, 0.03, 0.04]:
##    y = key_rate_EB(x, t, e)
##    y[y <  0] = 0
##    plt.plot(x,y)
#
#
#
##################### CHECK PM
#t = 0.1
#e = 0.005
#V = 1
#print("I_shared: ", I_shared_PM(V, t, e, 2))
#print("I_shared2:", I_shared_PM2(V, t, e, 2))
#
#
#print("Key rate PM:", key_rate_PM(V, t, e))
#print("Key rate PM:", key_rate_PM(V, t, e, 2))
#
######## PAPER figures validation
# plt.figure()
# x = np.linspace(0.001, 50, 5000)
# t = .1
# for e in [0.0, 0.001, 0.002, 0.003, 0.004]:
#    y = key_rate_PM(x, t, 2*e)
#    y[y <  0] = 0
#    plt.plot(x,y)
#
#
# plt.figure()
# t = .1
# e = 0.005
# x = np.linspace(0.001, 30, 5000)
# for b in [1, 0.98, 0.96, 0.94, 0.92, 0.90]:
#    y = key_rate_PM_beta(x, t, e, b)
#    y[y <  0] = 0
#    plt.plot(x,y)
#




#
# plt.figure()
# L = np.linspace(0, 100)
# t = L *.2
# t = 10**(-t/10)
#
# for e in [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]:
#     V_opt, k_opt = opt_key_rate(t, e, 2)
#     k_opt[k_opt < 0] = 0
#     plt.plot(L, k_opt)


# plt.figure()
# x = np.linspace(0.001, 5000, 5000)
# y = key_rate_PM(x, 1, 0, 1)
# plt.plot(x,y)
#
# y = key_rate_PM(x, 1, 0, 2)
# plt.plot(x,y)


# plt.figure()
# e = np.linspace(0, 0.1)
# for l in [5, 10, 20, 40]:
#     t = l * .2
#     t = 10**(-t/10)
#
#     V_opt, k_opt = opt_key_rate_e(t, 3*e, 2)
#     k_opt[k_opt < 0] = 0
#     plt.plot(e, k_opt)


#plt.figure()
#b = np.linspace(0.5, 1)
#t = 0.1
#for e in [0.001, 0.005, 0.01]:
#        V_opt, k_opt = opt_key_rate_b(t, e, b, 1)
#        k_opt[k_opt < 0] = 0
#        plt.plot(b, k_opt)
#
#plt.show()


# plt.figure()
# b = np.linspace(0.5, 1)
# e = 0.005
# for l in [5, 10, 20, 40]:
#         t = l * .4
#         t = 10**(-t/10)
#
#         V_opt, k_opt = opt_key_rate_b(t, e, b, 1)
#         k_opt[k_opt < 0] = 0
#         plt.plot(b, k_opt)
#
# plt.show()
