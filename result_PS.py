# -*- coding: utf-8 -*-
"""
Replicaing Photon subtraction results of Minjgian's paper

Created on Mon Apr  1 14:12:46 2019

@author: Eduardo Villasenor
"""

import cv_system as cv
import measurements
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.io


############################################ CALCULATIONS

# Parameters
N = 20
mpn = 1.3
mpne = 0.001
f = 0.95
option = 'rps'

# Photon subtraction options
t = .9

# Scissors options
# Best  k=0.5 m_aux=0.6
k = .4
m_aux = .2



ps_theta = np.arccos(np.sqrt(t))
r = np.arcsinh(np.sqrt(mpn))
r_eve = np.arcsinh(np.sqrt(mpne))
r_aux = np.arcsinh(np.sqrt(m_aux))

## Initialize state
sys = cv.System(N, Nmodes=2, cm=False)
sys.apply_TMS(r, [0, 1])
# BAD TMSV
#sys.replace_current_state_w_bad_TMSV(mean_photon_number)


# Transmitter Photon subtraction
if option == 'tps':
    p_success = sys.apply_photon_subtraction(t, 1)
    print("P SUCCESS:", p_success)

elif option == 'tct':
    p_success = sys.apply_photon_catalysis(1, t, 1)
    print("P SUCCESS:", p_success)

# No Photon subtraction
elif option == 'none':
    p_success = 1

# Transmitter Scissors
elif option == 'tsc':
#    p_success = sys.apply_scissors(k, r_aux, 1)
    p_success = sys.apply_scissors_options(k, r_aux, 1, 'c')
    print("P SUCCESS:", p_success)
#    print(sys.state)

# Transmitter Scissors Exact
elif option == 'tsc_e':
    p_success = sys.apply_scissors_exact(k, 1)
    print("P SUCCESS:", p_success)
#    print(sys.state)


# Evesdropper collective attack
sys.add_TMSV(r_eve)
# BAD TMSV
#sys.add_bad_TMSV(e_mpn)

sys.set_quadratures_basis()


# Save current state of the system
sys.save_state()

key_rates = []
tes = np.logspace(-2, 0, base=10, num=20)
#tes = np.linspace(.9, 1, 10)
#tes = [1.]

for te in tes:
    sys.load_state()

    theta = np.arccos(np.sqrt(te))
    sys.apply_BS(theta, [1, 2])

    # Receiver Photon subtraction
    if option == 'rps':
        p_success = sys.apply_photon_subtraction(t, 1)
        print("P SUCCESS:", p_success)
        
    elif option == 'rct':
        p_success = sys.apply_photon_catalysis(1, t, 1)
        print("P SUCCESS:", p_success)

    # Receiver Scissors
    elif option == 'rsc':
        p_success = sys.apply_scissors(k, r_aux, 1)
        print("P SUCCESS:", p_success)

    # Receiver Scissors Exact
    elif option == 'rsc_e':
        p_success = sys.apply_scissors_exact(k, 1)
        print("P SUCCESS:", p_success)

#    print(sys.cm)
#    print(sys.get_full_CM())
    kr = measurements.key_rate(sys, f, p_success)
#    kr = measurements.key_rate_compare(sys, f, p_success, mpn, mpne, te)
    key_rates += [kr]
    print("--->", te, kr)


# Save the resuls
filename = "data/result_PS_" + option
key_rates = np.array(key_rates)
#print(key_rates)
np.save(filename, key_rates)

filename_ind = "data/indeces_PS_" + option
np.save(filename_ind, tes)


############################################ PLOT

key_rates = []
indeces = []
for ext in ['none', 'tps', 'rps', 'tsc_e', 'rsc_e', 'tct', 'rct']:
    f_name = "data/result_PS_" + ext
    k_rate = np.load(f_name + ".npy")
    key_rates += [k_rate]

    i_name = "data/indeces_PS_" + ext
    inds = np.load(i_name + '.npy')
    indeces += [inds]


fig1, ax = plt.subplots()
lines_types = ['k*-', 'b*-', 'r*-', 'y*-', 'm*-', 'c*-', 'g*-']
for i in range(len(key_rates)):
    ax.plot(indeces[i], key_rates[i], lines_types[i])


mg_ratesNO = scipy.io.loadmat('data/rate_no.mat')
mg_ratesR = scipy.io.loadmat('data/rate_r.mat')
mg_ratesT = scipy.io.loadmat('data/rate_t.mat')
mg_indNO = scipy.io.loadmat('data/ind_no.mat')
mg_indR = scipy.io.loadmat('data/ind_r.mat')
mg_indT = scipy.io.loadmat('data/ind_t.mat')

mg_ratesNO = mg_ratesNO['RateNo'][0]
mg_ratesR = mg_ratesR['RateR'][0]
mg_ratesT = mg_ratesT['RateT'][0]
mg_indNO = mg_indNO['indd'][0]
mg_indR = mg_indR['indR'][0]
mg_indT = mg_indT['indT'][0]

ax.plot(mg_indNO, mg_ratesNO, 'k--')
ax.plot(mg_indR, mg_ratesR, 'r--')
ax.plot(mg_indT, mg_ratesT, 'b--')

ax.set_xlabel(r"$\eta$")
ax.set_ylabel("Key rate")
#ax.legend(["No PS", "Transmitter PS", "Receiver PS", "Transmitter SC", "Mg no PS", "Mg r-PS", "Mg t-PS"])


#key_rates1b = np.load(filename1 + "BAD.npy")
#key_rates2b = np.load(filename2 + "BAD.npy")
#key_rates3b = np.load(filename3 + "BAD.npy")
#ax.plot(tes, key_rates1b, 'ko')
#ax.plot(tes, key_rates2b, 'bo')
#ax.plot(tes, key_rates3b, 'ro')

ax.set_yscale('log')


locmaj = matplotlib.ticker.LogLocator(base=10,numticks=10)
ax.yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

#ax.minorticks_on()
plt.show()
