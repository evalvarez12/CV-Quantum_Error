import fidelities_DA as fda
import TMSV as tmsv
import coherent as co
import TMSV_PS as ps
import numpy as np
import matplotlib.pyplot as plt
import scipy.io



### Comparing fidelities with MGs code
tdb = np.linspace(0, 1, 15)
eps = 7e-3
sigma = 1
res = []
ps = []
for it in tdb:
    t =  it
    print(t)
    res += [tmsv.opt_fidelity_alphabet(t, eps, 1, 1, sigma)]


scipy.io.savemat('matlab/tmsv.mat', {'data':res})

#
#
#
#
#
#
#
# T = .6
# V = 6
# eps = 0.0
# alpha = 1 + 1j
# theta = 1E-6
# delta = .5
#
# # Ancilla definitions for DA
# Tda = 1
# eta = Tda
# if T == 1:
#     nth = 0
# else:
#     nth = eps/(2*(1 - T))
# r = np.arccosh(V)/2
# tau = -np.log(T)
# g = 1
#
# print('Compare TMSV')
# print('Direct:', co.fidelity(T, eps, alpha))
# print('TMSV:', tmsv.fidelity(V, T, eps, eta, alpha))
# print('TMSV old:', tmsv.fidelity_old(V, T, eps, alpha))
#
# # print('TMSV DA:', fda.tmsv_eq(r, tau, g, Tda, np.real(alpha), np.imag(alpha), nth))
# print('----------------')
# print('PS:', ps.fidelity(V, theta, T, eps, alpha))
# # print('SB:', fda.squeezed_bell_eq(r, delta, tau, g, Tda, np.real(alpha), np.imag(alpha), nth))
#
# # print('SB DA:', fda.fidelity_(r, tau, g, Tda, np.real(alpha), np.imag(alpha)))
#
# # print('PS2:', ps2.fidelity(V, theta, T, eps, alpha=1))
# print('--------- Optimized')
# # print('--> opt TMSV:', tmsv.opt_fidelity(T, eps, alpha))
# # print('--> opt PS:', ps.opt_fidelity(T, eps, alpha))
# # print('--> opt PS V', ps.opt_fidelity_V(3.99, T, eps, alpha))
# # print('--> opt TMSV DA:', fda.get_opt_f('tmsv', tau, g, Tda, alpha, nth))
# # print('--> opt SB:', fda.get_opt_f('squeezed_bell', tau, g, Tda, alpha, nth))
# print('----------------')
#
#
# ###################### PLOT
# plt.rcParams["font.family"] = "Times New Roman"
# plt.close('all')
# # plt.figure()
# # T = np.linspace(1, 0)
# # plt.plot(T, co.fidelity(T, 0, alpha*0), 'o-')
# # plt.plot(T, co.fidelity(T, 0.05, alpha*0), 'o-')
#
# plot = True
# if plot:
#     plt.figure()
#     alpha = 1 + 1j
#     eps = 0.05
#
#     # Ancilla definitions for DA
#     Tda = 1
#     g = 1/Tda
#     eta = Tda
#
#     F_tm = []
#     F_ps = []
#
#     F_sbDA = []
#     F_scDA = []
#     F_tmDA = []
#
#     trans = np.linspace(0.1,.99, 20)
#
#     for ti in trans:
#         print('--------------', ti)
#         print('------------> TMSV')
#         F_tm += [tmsv.opt_fidelity(ti, eps, eta, alpha)]
#         F_ps += [ps.opt_fidelity(ti, eps, alpha)]
#
#         tau = -np.log(ti)
#         if ti == 1:
#             nth = 0
#         else:
#             nth = eps/(2*(1 - ti))
#
#         # print('------------> TMSV DA')
#         # F_tmDA += [fda.get_opt_f('tmsv', tau, g, Tda, alpha, nth)]
#         print('------------> SB DA')
#         F_sbDA += [fda.get_opt_f('squeezed_bell', tau, g, Tda, alpha, nth)]
#         # F_sc += [fidelity.get_opt_f_r('squeezed_cat', ri, t, g, T, B)]
#
#
#     plt.plot(trans, F_tm, 'k--', label='TMSV')
#     plt.plot(trans, F_ps, 'v', label='PS')
#
#
#     # plt.plot(trans, F_tmDA, 's' ,label='TMSV_DA')
#     plt.plot(trans, F_sbDA, 'o', label='SB')
#
#
#     plt.grid()
#     plt.legend()
#     plt.xlabel('Transmissivity')
#     plt.ylabel('Fidelity')
#
#     # plt.ylim((.45,1))
#
#
# plot = True
# if plot:
#     plt.figure()
#     alpha = 1 + 1j
#     eps = 0.05
#
#     # Ancilla definitions for DA
#     Tda = 1
#     g = 1/Tda
#     eta = Tda
#
#     V = np.linspace(2,7)
#
#     # Fix Transmissivity
#     T = [0.5]
#
#     for Ti in T:
#         F_tm = []
#         F_ps = []
#
#         F_sbDA = []
#         F_scDA = []
#         F_tmDA = []
#
#         # Alternative parameters defs
#         tau = -np.log(Ti)
#         nth = eps/(2*(1- Ti))
#
#         for Vi in V:
#             # My functions
#             # print('--------------', ti)
#             # print('------------> TMSV')
#             F_tm += [tmsv.fidelity(Vi, Ti, eps, eta, alpha)]
#             F_ps += [ps.opt_fidelity_V(Vi, Ti, eps, alpha)]
#
#
#             # DA functions
#             r = np.arccosh(Vi)/2
#             # print('------------> TMSV DA')
#             # F_tmDA += [fda.get_opt_f('tmsv', tau, g, Tda, alpha, nth)]
#             # print('------------> SB DA')
#             F_sbDA += [fda.get_opt_f_r('squeezed_bell', r, tau, g, Tda, alpha, nth)]
#             # F_sc += [fidelity.get_opt_f_r('squeezed_cat', ri, t, g, T, B)]
#
#
#         plt.plot(V, F_tm, 'k--', label='TMSV')
#         plt.plot(V, np.array(F_ps), 'v-', label='PS')
#
#         # plt.plot(V, F_tmDA, 's' ,label='TMSV_DA')
#         plt.plot(V, F_sbDA, 'o-', label='SB')
#
#         plt.grid()
#         plt.legend()
#         plt.xlabel('V')
#         plt.ylabel('Fidelity')

        # plt.ylim((.45,1))




plt.show()


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
