# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 14:24:13 2019

@author: z5239621
"""

import TMSV as tmsv
import coherent as co
import scipy.io
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

# Global extra loss parameters
absorption = lambda z : np.exp(-0.7/np.cos(np.deg2rad(z)))
#absorption = lambda z : 1
#eta_db = 1
#eta = 10**(-eta_db/10)
eta = 0.9
t_ext_db = 1
t_ext = 10**(-t_ext_db/10)


############ UPLINK & DOWNLINK CHANNELS PROPERTIES PLOT


#
# plt.figure()
# #for z in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]:
# r = 1
#
# colors_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
# colors = {0: colors_cycle[0], 15:colors_cycle[1], 30:colors_cycle[2], 45:colors_cycle[3], 60:colors_cycle[4], 75:colors_cycle[5]}
#
# bins = np.arange(0.03, .450, .001)
# bins2 = np.arange(0.0, 1, .001)
#
# for z in [0, 15, 30, 45, 60]:
#     tag = "_z=" + str(z) + "_1024_10000"
#     data_file = "../../Laser_propagation/data/TELE_DOWN_I_r=" + str(r) + tag
#     downlink_data = scipy.io.loadmat(data_file)['res'].transpose()
#
#     data_file = "../../Laser_propagation/data/TELE_UP_I_r=" + str(r) + tag
#     uplink_data = scipy.io.loadmat(data_file)['res'].transpose()
#
#     print('z = ', z)
#     plt.hist(uplink_data * t_ext, bins=bins, histtype='step', label=str(z), density=True, linewidth=2, color=colors[z], linestyle=':')
#
# #    plt.hist((1-P)/P * db2, bins=bins2, histtype='step', label=str(z), density=True, linewidth=1.5, color=colors[z])


# for z in [0, 30, 75]:
#     tag = "_z=" + str(z) + "_1024_10000"
#     data_file = "l0d=0.1SH_DOWN_I_r=" + str(r)  + tag
#     I = scipy.io.loadmat(data_file)['res']
#     plt.hist(I * t_ext, bins=bins, histtype='step', label=str(z), density=True, linewidth=1.5, linestyle=':', color=colors[z])


#    plt.hist(P0*I, bins=bins)
    #plt.hist(I, bins=bins)
plt.xlabel('$T$')
plt.ylabel('PDF')
plt.legend()


scint_up_dic = {}
scint_down_dic = {}

colors = {'uplink':'b', 'downlink':'g'}
zs = [0, 5, 10, 15, 20 , 25, 30, 35, 40, 45, 50, 55, 60]
r = 1
fig = plt.figure()
ax = fig.add_subplot(111)
# ax2 = ax.twinx()
uplink_avg = []
uplink_std = []
downlink_avg = []
downlink_std = []

downlink_eff = []
uplink_eff = []
downlink_eff_noise = []
uplink_eff_noise = []

for z in zs:

    tag = "_z=" + str(z) + "_1024_10000"
    data_file = "../../Laser_propagation/data/TELE_final2_DOWN_I_r=" + str(r) + tag
    downlink_data = scipy.io.loadmat(data_file)['res'].transpose()

    data_file = "../../Laser_propagation/data/TELE_final2_UP_I_r=" + str(r) + tag
#    data_file = "../../Laser_propagation/data/TELE__L0inf_UP_I_r=" + str(r) + tag
    uplink_data = scipy.io.loadmat(data_file)['res'].transpose()


    T_down = downlink_data * t_ext * absorption(z)
    T_up = uplink_data * t_ext * absorption(z)


    # Define effective parameters
    T_eff_down = np.average(np.sqrt(T_down))**2
    scint_down = np.average(T_down**2)/(np.average(T_down)**2) - 1 + 0.007
    eps_eff_down = np.var(np.sqrt(T_down))/T_eff_down + np.average(T_down)/T_eff_down * scint_down

    T_eff_up = np.average(np.sqrt(T_up))**2
    scint_up = np.average(T_up**2)/(np.average(T_up)**2) - 1 + 0.007
    eps_eff_up = np.var(np.sqrt(T_up))/T_eff_up + np.average(T_up)/T_eff_up * scint_up

    # Convert to dB
    T_eff_up = -10 * np.log10(T_eff_up)
    T_eff_down = -10 * np.log10(T_eff_down)

    T_down_db = -10 * np.log10(T_down)
    T_up_db = -10 * np.log10(T_up)
#    print('z=', z, 'avg', np.round(-10 * np.log10(np.average(T_up)), 3), 'std', np.round(-10 * np.log10(np.std(T_up)),3))
#    print('z=', z, 'avg', np.average(T_up_db), 'std', np.std(T_up_db))
#    print(z, np.round(-10 * np.log10(np.average(T_up)), 3))
#    print(z, np.round(-10 * np.log10(np.average(T_down)), 3))


#    print('z=', z, 'avg', -10 * np.log10(np.average(T_down)), 'std', -10 * np.log10(np.std(T_down)))
#    print('z=', z, 'avg', np.average(T_down_db), 'std', np.std(T_down_db))

    downlink_avg += [np.average(T_down_db)]
    downlink_std += [np.std(T_down_db)]

    downlink_eff += [T_eff_down]
    downlink_eff_noise += [eps_eff_down]

    uplink_avg += [np.average(T_up_db)]
    uplink_std += [np.std(T_up_db)]

    uplink_eff += [T_eff_up]
    uplink_eff_noise += [eps_eff_up]


#    print('z=', z, 'avg', np.round(uplink_avg[-1],3), 'std', np.round(uplink_std[-1],3))
#       print('z=', z, 'avg', np.round(downlink_avg[-1],3), 'std', np.round(downlink_std[-1],3))

    scint_down_dic.update( {z : np.average(T_down**2)/(np.average(T_down)**2) - 1 + 0.007})
    scint_up_dic.update({z : np.average(T_up**2)/(np.average(T_up)**2) - 1 + 0.007})







############################## UPLINK & TELEPORTATION FIDELITIES PLOT
fig = plt.figure()
ax = fig.add_subplot(111)

zs = [0, 5, 10, 15, 20 , 25, 30, 35, 40, 45, 50, 55, 60]
r = 1

sigmaT = [2, 4, 10, 25]
sigmaD = [2, 4, 6]


uplink_avg = []
uplink_std = []
downlink_avg = []
downlink_std = []

f_up1_avg = []
f_up1_std = []
f_up2_avg = []
f_up2_std = []
f_up3_avg = []
f_up3_std = []
f_up4_avg = []
f_tmsv_avg = []
f_tmsv_std = []

f_tmsv_avg1 = []
f_tmsv_avg2 = []
f_tmsv_avg3 = []
f_tmsv_avg4 = []


for z in zs:

    tag = "_z=" + str(z) + "_1024_10000"
    data_file = "../../Laser_propagation/data/TELE_final2_DOWN_I_r=" + str(r) + tag
    downlink_data = scipy.io.loadmat(data_file)['res'].transpose()

    data_file = "../../Laser_propagation/data/TELE_final2_UP_I_r=" + str(r) + tag
    uplink_data = scipy.io.loadmat(data_file)['res'].transpose()


    T_down = downlink_data * t_ext * absorption(z)
    T_up = uplink_data * t_ext * absorption(z)


    # Define effective parameters
    T_down_avg = np.average(T_down)
    eps_down =  scint_down_dic[z]

    # print('T_f DOWN',z,T_down_avg)
#    print('e_f DOWN',z,eps_down)
    T_up_avg = np.average(T_up)
    eps_up = scint_up_dic[z]
#    print('T_f UP',z,T_up)
#    print('e_f UP',z,eps_up)
   # TMSV
#       V_opt, F_opt = tmsv.opt_values(T_eff_down, eps_eff_down, eta, alpha=sigma[2])
#       print("z =", z, V_opt, F_opt)
#

    pars_t_1 = tmsv.opt_avg_fidelity_vareps_getoptpars(T_down_avg, eps_down, eta, sigmaT[0])
    pars_t_2 = tmsv.opt_avg_fidelity_vareps_getoptpars(T_down_avg, eps_down, eta, sigmaT[1])
    pars_t_3 = tmsv.opt_avg_fidelity_vareps_getoptpars(T_down_avg, eps_down, eta, sigmaT[2])
    pars_t_4 = tmsv.opt_avg_fidelity_vareps_getoptpars(T_down_avg, eps_down, eta, sigmaT[3])



    f_tmsv_avg1 += [np.average(tmsv.fidelity_alphabet_pars_vareps(pars_t_1, T_down, eps_down, eta, sigmaT[0]))]
    f_tmsv_avg2 += [np.average(tmsv.fidelity_alphabet_pars_vareps(pars_t_2, T_down, eps_down, eta, sigmaT[1]))]
    f_tmsv_avg3 += [np.average(tmsv.fidelity_alphabet_pars_vareps(pars_t_3, T_down, eps_down, eta, sigmaT[2]))]
    f_tmsv_avg4 += [np.average(tmsv.fidelity_alphabet_pars_vareps(pars_t_4, T_down, eps_down, eta, sigmaT[3]))]


    eps = eps_up * sigmaD[0]
    f_up1_avg += [np.average(co.fidelity_alphabet(T_up, eps, sigmaD[0]))]


    eps = eps_up * sigmaD[1]
    f_up2_avg += [np.average(co.fidelity_alphabet(T_up, eps, sigmaD[1]))]

    eps = eps_up * sigmaD[2]
    f_up3_avg += [np.average(co.fidelity_alphabet(T_up, eps, sigmaD[2]))]

#       eps = eps_eff_up * sigma[3]
#       f_dt4 = co.fidelity_alphabet(T_eff_up, eps, sigmaD[3])
#       f_up4_avg += [f_dt4]

#   ax.errorbar(zs, f_tmsv_avg, f_tmsv_std, label=r'Teleportation', linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5)
#   ax.errorbar(zs, f_up1_avg, f_up1_std, label=r'Uplink $\sigma= $' + str(np.round(np.abs(alpha[0]),2)), linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5)
#   ax.errorbar(zs, f_up2_avg, f_up2_std, label=r'Uplink $\sigma= $' + str(np.round(np.abs(alpha[1]),2)), linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5)
#   ax.errorbar(zs, f_up3_avg, f_up3_std, label=r'Uplink $\sigma= $' + str(np.round(np.abs(alpha[2]),2)), linestyle='--', marker='o', capsize=4, markersize=4, linewidth=1.5)

print(f_tmsv_avg4)

ax.plot(zs, f_tmsv_avg1, label=r'Upload $\sigma= $' + str(np.round(np.abs(sigmaT[0]),2)), linestyle='--', marker='o', markersize=4, linewidth=1.5)
ax.plot(zs, f_tmsv_avg2, label=r'Upload $\sigma= $' + str(np.round(np.abs(sigmaT[1]),2)), linestyle='--', marker='o', markersize=4, linewidth=1.5)
ax.plot(zs, f_tmsv_avg3, label=r'Upload $\sigma= $' + str(np.round(np.abs(sigmaT[2]),2)), linestyle='--', marker='o', markersize=4, linewidth=1.5)
ax.plot(zs, f_tmsv_avg4, label=r'Upload $\sigma= $' + str(np.round(np.abs(sigmaT[3]),2)), linestyle='--', marker='o', markersize=4, linewidth=1.5)


ax.plot(zs, f_up1_avg, label=r'Uplink $\sigma= $' + str(np.round(np.abs(sigmaD[0]),2)), linestyle='--', marker='*', markersize=4, linewidth=1.5)
ax.plot(zs, f_up2_avg, label=r'Uplink $\sigma= $' + str(np.round(np.abs(sigmaD[1]),2)), linestyle='--', marker='*', markersize=4, linewidth=1.5)
ax.plot(zs, f_up3_avg, label=r'Uplink $\sigma= $' + str(np.round(np.abs(sigmaD[2]),2)), linestyle='--', marker='*', markersize=4, linewidth=1.5)
#ax.plot(zs, f_up4_avg, label=r'Uplink $\sigma= $' + str(np.round(np.abs(sigma[3]),2)), linestyle='--', marker='o', markersize=4, linewidth=1.5)




#plt.plot(zs, np.zeros_like(zs), ls='--', c='k')

clasical = np.ones_like(f_tmsv_avg1) * .5
plt.plot(zs, clasical, 'r--')


ax.set_xlabel(r'$\zeta$ (deg)')
ax.set_ylabel(r'$\langle \mathcal{F} \rangle $')
# plt.ylim([-0.005, 0.2])
#ax.set_ylim(0.495, 0.65)
ax.set_xlim(-0.5, 60.5)
ax.grid()
ax.legend()

### To MATLAB
data = [zs, f_tmsv_avg1, f_tmsv_avg2, f_tmsv_avg3, f_tmsv_avg4, f_up1_avg, f_up2_avg ,f_up3_avg]
scipy.io.savemat('matlab/fidelities.mat', {'data':data})

plt.show()
