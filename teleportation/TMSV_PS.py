import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt


def p_success(V, theta, T, eps, alpha):
    # Define elements of covariance matrix
    a = V
    c = np.sqrt(T * (V**2 - 1))
    b = T * (V - 1) + 1 + eps

    # Beam splitter properties
    t = np.cos(theta)
    r = np.sin(theta)

    A = (b*r**2 + t**2 + 1)/2
    return 1/A - 1/A**2


def fidelity(V, theta, T, eps, alpha):
    gp = 1
    gx = 1

    # Define elements of covariance matrix
    a = V
    c = np.sqrt(T * (V**2 - 1))
    b = T * (V - 1) + 1 + eps

    # print('a:', a)
    # print('b:', b)
    # print('c:', c)

    # Beam splitter properties
    t = np.cos(theta)
    r = np.sin(theta)

    # print('t:', t)
    # print('r:', r)

    # First definitions
    A = (b*r**2 + t**2 + 1)/2

    A2u = (a*gp**2 + b*t**2 + r**2)/2 - gp*c*t - r**2*((b-1)*t - gp*c)**2/(4*A) + (gp+1)/2
    A2v = (a*gx**2 + b*t**2 + r**2)/2 - gx*c*t - r**2*((b-1)*t - gx*c)**2/(4*A) + (gx+1)/2

    B2u = -2*np.imag(alpha)*(1-gp)
    B2v = 2*np.real(alpha)*(1-gx)

    D1u = r**2*((b-1)*t - c*gp)**2 / (4 * A**3)
    D1v = r**2*((b-1)*t - c*gx)**2 / (4 * A**3)

    AA = 1/A - 1/A**2

    E = np.exp(-B2u**2/(4*A2u) - B2v**2/(4*A2v))

    #######
    # print('A:', A)
    # print('A2u:', A2u)
    # print('A2v:', A2v)
    # print('B2u:', B2u)
    # print('B2v:', B2v)
    # print('D1u:', D1u)
    # print('D1v:', D1v)
    # print('E:', E)
    # print('p:', AA)
    ########

    res = (AA/np.sqrt(A2u * A2v) - D1u * (2 * A2u - B2u**2)/(4 * np.sqrt(A2u**5 * A2v)) - D1v * (2 * A2v - B2v**2)/(4* np.sqrt(A2v**5 * A2u)) ) * E/(AA)
    return res

def fidelity_pars(pars, T, eps, alpha):
    V, theta = pars
    return fidelity(V, theta, T, eps, alpha)

def opt_fidelity(T, eps, alpha):
    F = lambda pars : 1 - fidelity_pars(pars, T, eps, alpha)
    initial_guess = [2, 2]
    cons=({'type': 'ineq',
       'fun': lambda x: x},
      {'type': 'ineq',
       'fun': lambda x: 2 * np.pi  - np.abs(x)})

    res = op.minimize(F, initial_guess, constraints=cons)
    if not res['success']:
        print(res)
        # raise AssertionError('Failure in optimization')
    return fidelity_pars(res['x'], T, eps, alpha)

# print(opt_fidelity(1, 0.001, 1))
# alpha = 1
# V = 6
# theta = np.pi/4
#
# eps = 0.001
# T = np.linspace(0, 1)
# plt.plot(T, fidelity(V, theta, T, eps, alpha))
# plt.show()


# print(fidelity(5, np.pi/2, 1, 0, 1))
