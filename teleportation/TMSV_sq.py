import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import math

def fidelity(V, T, eps, g, eta, s):
    g = 1/eta
    A = (V + T*(V-1) + 1 + eps)*(g*eta)**2 + g**2*(1-eta**2) - 2 *g*eta*np.sqrt(T*(V**2-1))
    C = (g*eta)**2 +1

    F = 2/np.sqrt((A+C*(np.cosh(s) + np.sinh(s))**2)*(A+C*(np.cosh(s) - np.sinh(s))**2))
    # if not math.isnan(F):
    #     print('Vg',V, g)
    #     print('FA:', F,A)

    return F

def fidelity_pars(pars, T, eps, eta, s):
    V, g = pars
#    g = 1/eta
    return fidelity(V, T, eps, g, eta, s)


def opt_fidelity(T, eps, eta, s):
    F = lambda P : 1 - fidelity_pars(P, T, eps, eta, s)
    initial_guess = [1.2, 1.1]
    cons=({'type': 'ineq',
        'fun': lambda x: x[0] - 1.001},
          {'type': 'ineq',
        'fun': lambda x: x[1]})
    # res = op.minimize(F, initial_guess, constraints=cons)
    res = op.minimize(F, initial_guess, constraints = cons)
    print(res)
    # if not res['success']:
        # raise AssertionError('Failure in optimization')
    return fidelity_pars(res['x'], T, eps, eta, s)


def opt_fidelity_params(T, eps, eta, s):
    F = lambda P : 1 - fidelity_pars(P, T, eps, eta, s)
    initial_guess = [1.2, 1.1]
    cons=({'type': 'ineq',
        'fun': lambda x: x[0] - 1.001},
          {'type': 'ineq',
        'fun': lambda x: x[1]})
    # res = op.minimize(F, initial_guess, constraints=cons)
    res = op.minimize(F, initial_guess, constraints = cons)
    print(res)
    # if not res['success']:
        # raise AssertionError('Failure in optimization')
#    return fidelity_pars(res['x'], T, eps, eta, s)
    return res['x']


def opt_fidelity_avg(T, eps, eta, s_mean, sigma):
    s = np.random.normal(s_mean, np.sqrt(sigma), 20000)
    s = np.abs(s)
    opt_params = opt_fidelity_params(T, eps, eta, s_mean)
    V_opt, g_opt = opt_params
    
#    print(opt_params, s)
#    print(V_opt, T, eps, eta, g_opt, s_mean)
#    print('F:', fidelity(V_opt, T, eps, g_opt, eta, s_mean))
#    print('F:', fidelity_pars(opt_params, T, eps, eta, s_mean))
    F = fidelity(V_opt, T, eps, g_opt, eta, s)
#    print(F)
    return np.average(F)