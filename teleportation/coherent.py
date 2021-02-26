import numpy as np
import scipy.optimize as op

def fidelity(T, eps, alpha):
    E = np.exp(-2 * (1 - np.sqrt(T))**2 * np.abs(alpha)**2 / (eps + 2))
    return  2 * E / (eps + 2)


def fidelity_amp(T, eps, alpha, g):
    eps = eps * g
    E = np.exp(-2 * (1 - g*np.sqrt(T))**2 * np.abs(alpha)**2 / (eps + 2))
    return  2 * E / (eps + 2)

def fidelity_alphabet(T, eps, sigma):
    F =  2 / (2*sigma*(1-np.sqrt(T))**2 + 2 + eps)
    return F

def fidelity_alphabet_amp(T, eps, sigma, g):
    eps = eps * g
    F =  2 / (2*sigma*(1-g*np.sqrt(T))**2 + 2 + eps)
    return F



def opt_fidelity_alphabet(T, eps, sigma):
    F = lambda g : 1 - fidelity_alphabet_amp(T, eps, sigma, g)
    initial_guess = 1
    cons=({'type': 'ineq',
       'fun': lambda x: x})
    res = op.minimize(F, initial_guess, constraints=cons)
#    res = op.minimize(F, initial_guess)
    print('res:', res)
    # if not res['success']:
        # raise AssertionError('Failure in optimization')
#    print('opt V:', np.round(res['x'],3))
    return fidelity_alphabet_amp(T, eps, sigma, res['x'])