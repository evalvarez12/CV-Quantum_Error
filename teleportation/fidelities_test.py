import fidelities_DA as fda
import TMSV as tmsv
import TMSV_PS2 as ps2
import TMSV_PS as ps
import numpy as np
import matplotlib.pyplot as plt

# Rob comments:

T = np.random.rand()
V = 1
eps = np.random.rand()
alpha = 1 + 1j
theta = .1

# Ancilla definitions for DA
Tda = 1
if T == 1:
    nth = 0
else:
    nth = eps/(2*(1 - T))
r = np.arccosh(V)/2
tau = -np.log(T)
g = 1

a = V
b = T*(V-1) + 1 + eps
c = np.sqrt(T * (V**2 - 1))
ref = 2 / (a + b - 2*c + 2)
Delta = fda.Delta(r, tau, g, Tda, nth)
print('Ref:', ref)
print('Delta:', Delta)
print('Compare TMSV')
print('TMSV:', tmsv.fidelity(V, T, eps, alpha))
print('TMSV DA:', fda.tmsv_eq(r, tau, g, Tda, np.real(alpha), np.imag(alpha), nth))
print('----------------')
print('PS:', ps.fidelity(V, theta, T, eps, alpha))
# print('SB DA:', fda.fidelity_(r, tau, g, Tda, np.real(alpha), np.imag(alpha)))

# print('PS2:', ps2.fidelity(V, theta, T, eps, alpha=1))
# print('--------- Optimized')
# print('opt TMSV:', tmsv.opt_fidelity(T, eps, alpha))
# print('opt TMSV DA:', fda.get_opt_f('tmsv', tau, g, Tda, alpha, nth))
# print('----------------')


###################### PLOT
# t = 0.05*0
# nth = 0.1
# T =  1
# B = 1 + 1j
# g = 1/T
# eps = nth
#
# alpha = 1 + 1j
#
# F_tm = []
# F_ps = []
#
# F_sbDA = []
# F_scDA = []
# F_tmDA = []
#
# trans = np.linspace(0.1,1)
#
# for ti in trans:
#     F_tm += [tmsv.opt_fidelity(ti, eps, alpha)]
#     F_ps += [ps.opt_fidelity(ti, eps, alpha)]
#
#
#
# tda = -np.log(trans)
# for ti in tda:
#     F_sbDA += [fidelity.get_opt_f('squeezed_bell', ti, g, T, B)]
#     # F_sc += [fidelity.get_opt_f_r('squeezed_cat', ri, t, g, T, B)]
#     F_tmDA += [fidelity.get_opt_f('tmsv', ti, g, T, B)]
#
#
# plt.plot(trans, F_tm, label='TMSV')
# plt.plot(trans, F_ps, label='PS')
#
# plt.plot(trans, F_tmDA, 's' ,label='TMSV_DA')
# plt.plot(trans, F_sbDA, 'o', label='SB')
# plt.rcParams["font.family"] = "Times New Roman"
#
# plt.legend()
# plt.xlabel('Transmissivity')
# plt.ylabel('Fidelity')
#
# plt.show()


###############################################################################

# t = 0
# nth = 0
# r = 1
# T =  1
# B = 1 + 1j
# g = 1

# res1 = fidelity.opt_tmsv(t, g, T, B)
# res2 = fidelity.opt_squeezed_cat(t, g, T, B)
# res3 = fidelity.opt_squeezed_bell(t, g, T, B)
#
# res2r = fidelity.opt_squeezed_cat_r(r, t, g, T, B)
# res3r = fidelity.opt_squeezed_bell_r(r, t, g, T, B)



# r = np.random.rand() * 10
# d = np.random.rand() * np.pi
# gm = np.random.rand() * 5
#
# r = 6
# d = 2.5
# gm =0.5
#
# print(r,d,gm)
# print(fidelity.squeezed_cat_eq(r, d, gm, .0, 1, 1, 0, 1))
# print(fidelity.squeezed_cat_eq2(r, d, gm, .0, 1, 1, 0, 1))
# print(fidelity.squeezed_cat_red(r, d, gm))
#
# print("---------------")
# a1 = np.random.rand()
# a2 = np.random.rand()
# a3 = np.random.rand()
# a4 = np.random.rand()
# a5 = np.random.rand()
# print(fidelity.squeezed_cat_eq(r, d, gm, a1, a2, a3, a4, a5))
# print(fidelity.squeezed_cat_eq2(r, d, gm, a1, a2, a3, a4, a5))
# print(fidelity.squeezed_cat_red(r, d, gm))
# print((1/((np.exp(-2) + 1) * (1 + np.exp(-1)))) * (1 + np.exp(-1 + (1/(4 * (np.exp(-2) + 1))) * (np.exp(-2) - np.exp(2)))))



##########################

# V = np.linspace(1, 20)
# T = [1, .8, .6]
#
#
# alpha = 1
# eps = 0.
#
# for t in T:
#     plt.plot(V, tmsv.fidelity(V, t, eps, alpha))
#     plt.plot(V, ps.fidelity(V, 0, t, eps, alpha))


# print(fidelity(5, 1, 0, 1))



# r = np.arccosh(V)/2
# tau = -np.log(T)
# T =  1
# B = 1 + 1j
# g = 1/T
#
# F_tmsv = []
# F_sc = []
# for t in tau:
#     plt.plot(V, fidelity.tmsv_eq(r, t, g, T, np.real(B), np.imag(B)), 'o')
#
