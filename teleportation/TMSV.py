import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt


def fidelity(V, T, eps, alpha):
    gp = 1
    gx = gp

    # Define elements of covariance matrix
    a = V
    c = np.sqrt(T * (V**2 - 1))
    b = T * (V - 1) + 1 + eps
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


def fidelity_pars(pars, T, eps, alpha):
    V = pars
    return fidelity(V, T, eps, alpha)

def opt_fidelity(T, eps, alpha):
    F = lambda pars : 1 - fidelity_pars(pars, T, eps, alpha)
    initial_guess = 1
    cons=({'type': 'ineq',
       'fun': lambda x: x})
    # res = op.minimize(F, initial_guess, constraints=cons)
    res = op.minimize(F, initial_guess)
    print(res)
    # if not res['success']:
        # raise AssertionError('Failure in optimization')
    return fidelity_pars(res['x'], T, eps, alpha)

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
