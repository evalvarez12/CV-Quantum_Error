import fidelities_DA as fda
import TMSV as tmsv
import TMSV_PS2 as ps2
import coherent as co
import TMSV_PS as ps
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.close('all')


#p(x) = \frac{1}{\sqrt{ 2 \pi \sigma^2 }} e^{ - \frac{ (x - \mu)^2 } {2 \sigma^2} },
# Distribution properties
t_avg = .6
sigma = .05
n = 20000

# Channel
eps = 0.05
T = np.random.normal(t_avg, sigma, n)
alpha = 1 + 1j

# TMSV
V_opt, F_avg = tmsv.opt_values(t_avg, eps, alpha)
Ftmsv = tmsv.fidelity(V_opt, T, eps, alpha)

# plt.plot(T, "k.")
plt.figure()
plt.plot(T, Ftmsv, "k.")
plt.plot(np.linspace(min(T), max(T), n), F_avg * np.ones_like(Ftmsv), 'r--')


# TMSV - fully optimised
T = np.random.normal(t_avg, sigma, n)
Ftmsv_fo = []
for ti in T:
    Ftmsv_fo += [tmsv.opt_fidelity(ti, eps, alpha)]
# plt.plot(T, "k.")
plt.plot(T, Ftmsv_fo, "k.")
plt.plot(np.linspace(min(T), max(T), n), np.average(Ftmsv_fo) * np.ones_like(Ftmsv), 'r--')


# SB
T = np.random.normal(t_avg, sigma, n)

Tda = 1
tau = -np.log(T)
nth = eps/(2*(1 - T))
tau_avg = -np.log(t_avg)
nth_avg = eps/(2*(1 - t_avg))
g = 1


opt_vals, F_avg = fda.get_opt_values('squeezed_bell', tau_avg, g, Tda, alpha, nth_avg)
r_opt, d_opt = opt_vals
Fsb = []
for i in range(len(T)):
    taui = tau[i]
    nthi = nth[i]
    Fsb += [fda.squeezed_bell_eq(r_opt, d_opt, taui, g, Tda, np.real(alpha), np.imag(alpha), nthi)]

# plt.plot(T, "k.")
plt.plot(T, Fsb, "c.")
plt.plot(np.linspace(min(T), max(T), n), F_avg * np.ones_like(Fsb), 'y--')


plt.figure()
bins = np.linspace(.5, 1, 50)
Ht, Bt = np.histogram(T, bins, density=True)
Htm, Btm = np.histogram(Ftmsv, bins, density=True)
Hsb, Bsb = np.histogram(Fsb, bins, density=True)
Htmfo, Btmfo = np.histogram(Ftmsv_fo, bins, density=True)


# plt.plot(Bt[1:], Ht, 'k-')
plt.plot(Btm[1:], Htm, '-')
plt.plot(Bsb[1:], Hsb, '-')
plt.plot(Btmfo[1:], Htmfo, '-')

plt.show()
