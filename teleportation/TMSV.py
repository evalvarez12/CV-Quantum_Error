import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt


def fidelity(V, T, eps, eta, alpha):
    g = 1/eta
    g = g * eta
    R = 1 - eta**2

    a = V
    b = T * (V - 1) + 1 + eps
    c = np.sqrt(T * (V**2 - 1))

    A = (a + b * g**2 - 2*c*g)/2
    A = A + (g**2 + 1 + (R * g/eta)**2)/2

    Bu = -2 * np.imag(alpha) * (1 - g)
    Bv = 2 * np.real(alpha) * (1 - g)

    E = np.exp(-(Bu + Bv)**2/(4*A))

    return E/A


def fidelity_alphabet(V, T,eps, eta, g, sigma):
    g = g * eta
    R = 1 - eta**2

    a = V
    b = T * (V - 1) + 1 + eps
    c = np.sqrt(T * (V**2 - 1))

    D = a + b * g**2 - 2*c*g + g**2 + 1 + (R * g/eta)**2

    F = 2 / (2*sigma * (1 - g)**2 + D)
    return F


def fidelity_old(V, T, eps, alpha):
    gp = 1
    gx = gp

    # Define elements of covariance matrix
    a = V
    b = T * (V - 1) + 1 + eps
    c = np.sqrt(T * (V**2 - 1))
    # print('------ V =', V, 'T =', T, 'eps =', eps)

    # print('a:', a)
    # print('b:', b)
    # print('c:', c)

    A1u = (a * gp**2 + b - 2*c*gp)/2
    A1v = (a * gx**2 + b - 2*c*gx)/2

    A2u = (gp**2 + 1)/2
    A2v = (gx**2 + 1)/2

    A3u = A1u + A2u
    A3v = A1v + A2v

    B1u = -2 * np.imag(alpha) * (-gp + 1)
    B1v = 2 * np.real(alpha) * (-gx + 1)

    E = np.exp(-B1u**2/(4*A3u) - B1v**2/(4*A3v))

    # Alternative method based in DellAnno paper
    # res = 2 / (V + T * (V-1) - 2 * np.sqrt(T*(V**2-1)) + 3 + 2 * eps)
    # return res

    #######
    # print('A1u:', A1u)
    # print('A1v:', A1v)
    # print('A2u:', A2u)
    # print('A2v:', A2v)
    # print('A3u:', A3u)
    # print('A3v:', A3v)
    # print('B1u:', B1u)
    # print('B1v:', B1v)
    # print('E:', E)
    ########

    res = E/np.sqrt(A3u * A3v)
    return res


def opt_fidelity(T, eps, eta, alpha):
    F = lambda V : 1 - fidelity(V, T, eps, eta, alpha)
    initial_guess = 1
    cons=({'type': 'ineq',
       'fun': lambda x: x})
    # res = op.minimize(F, initial_guess, constraints=cons)
    res = op.minimize(F, initial_guess)
    # print(res)
    # if not res['success']:
        # raise AssertionError('Failure in optimization')
    return fidelity(res['x'], T, eps, eta, alpha)


def opt_fidelity_alphabet(T, eps, eta, g, sigma):
    F = lambda V : 1 - fidelity_alphabet(V, T, eps, eta, g, sigma)
    initial_guess = 1
    cons=({'type': 'ineq',
       'fun': lambda x: x})
    # res = op.minimize(F, initial_guess, constraints=cons)
    res = op.minimize(F, initial_guess)
    # print(res)
    # if not res['success']:
        # raise AssertionError('Failure in optimization')
    print('opt V:', np.round(res['x'],3))
    return fidelity_alphabet(res['x'], T, eps, eta, g, sigma)


def opt_fidelity_alphabet_vareps(T, eps, eta, g, sigma):
    F = lambda V : 1 - fidelity_alphabet(V, T, V*eps, eta, g, sigma)
    initial_guess = 1
    cons=({'type': 'ineq',
       'fun': lambda x: x})
    # res = op.minimize(F, initial_guess, constraints=cons)
    res = op.minimize(F, initial_guess)
    # print(res)
    # if not res['success']:
        # raise AssertionError('Failure in optimization')
#    print('opt V:', np.round(res['x'],3))
    return fidelity_alphabet(res['x'], T, eps, eta, g, sigma)



def opt_fidelity_alphabet_vareps_gopt(T, eps, eta, sigma):
    F = lambda P : 1 - fidelity_alphabet_pars_vareps(P, T, eps, eta, sigma)
    initial_guess = [1, 1]
    cons=({'type': 'ineq',
       'fun': lambda x: x})
    # res = op.minimize(F, initial_guess, constraints=cons)
    res = op.minimize(F, initial_guess)
    # print(res)
    # if not res['success']:
        # raise AssertionError('Failure in optimization')
    # print('opt V:', np.round(res['x'],3))
    return fidelity_alphabet_pars(res['x'], T, eps, eta, sigma)


def fidelity_alphabet_pars_vareps(pars, T, eps, eta, sigma):
    V, g = pars
    return fidelity_alphabet(V, T, V*eps, eta, g, sigma)


def opt_fidelity_alphabet_gopt(T, eps, eta, sigma):
    F = lambda P : 1 - fidelity_alphabet_pars(P, T, eps, eta, sigma)
    initial_guess = [1, 1]
    cons=({'type': 'ineq',
       'fun': lambda x: x})
    # res = op.minimize(F, initial_guess, constraints=cons)
    res = op.minimize(F, initial_guess)
    # print(res)
    # if not res['success']:
        # raise AssertionError('Failure in optimization')
    # print('gopt opt V:', np.round(res['x'],3))
    return fidelity_alphabet_pars(res['x'], T, eps, eta, sigma)


def fidelity_alphabet_pars(pars, T, eps, eta, sigma):
    V, g = pars
    return fidelity_alphabet(V, T, eps, eta, g, sigma)


def opt_values(T, eps, eta, alpha):
    F = lambda V : 1 - fidelity(V, T, V * eps, eta, alpha)
    initial_guess = 1
    cons=({'type': 'ineq',
       'fun': lambda x: x})
    # res = op.minimize(F, initial_guess, constraints=cons)
    res = op.minimize(F, initial_guess)
    if not res['success']:
        print(res)
        # raise AssertionError('Failure in optimization')
    return res['x'], fidelity(res['x'], T, eps, eta, alpha)



# # print(opt_fidelity(1, 0.001, 1))
# alpha = 1
# # T = .9
#
# eps = 0.5
# for T in [1, .8, .6]:
#     V = np.linspace(1, 20, 200)
#     plt.plot(V, fidelity(V, T, eps, alpha))
# plt.show()
#
#
# # print(fidelity(5, 1, 0, 1))
