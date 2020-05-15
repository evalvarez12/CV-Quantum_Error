import fidelities as fds
import numpy as np
import matplotlib.pyplot as plt


fidelity = fds.Fidelity()


t = 0
nth = 0
r = 1
T =  1
B = 1 + 1j
g = 1

res1 = fidelity.opt_tmsv(t, g, T, B)
res2 = fidelity.opt_squeezed_cat(t, g, T, B)
res3 = fidelity.opt_squeezed_bell(t, g, T, B)


r = np.linspace(0, 3)
t = 0
nth = 0
T =  1
B = 1 + 1j
g = 1

F_sb = []
F_sc = []
F_tm = []

for ri in r:
    F_sb += [fidelity.get_opt_f_r('squeezed_bell', ri, t, g, T, B)]
    F_sc += [fidelity.get_opt_f_r('squeezed_cat', ri, t, g, T, B)]
    F_tm += [fidelity.get_opt_f_r('tmsv', ri, t, g, T, B)]
