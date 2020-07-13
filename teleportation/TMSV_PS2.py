import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt

# This is the old version one, here the collapse to the single photon CF includes
# the variable change of the BS during PS

def fidelity(V, theta, T, eps, alpha):
    gp = 1
    gx = 1

    # Define elements of covariance matrix
    a = V
    c = np.sqrt(T * (V**2 - 1))
    b = T * (V - 1) + 1 + eps

    print('a:', a)
    print('b:', b)
    print('c:', c)

    # Beam splitter properties
    t = np.cos(theta)
    r = np.sin(theta)

    print('t:', t)
    print('r:', r)

    # First definitions
    A = b * r**2/2 + t**2

    A2u = gp**2 * a/2 + b*t**2/2 + r**2 - c*t*gp + ((b-2)*t*r - c*gp*r)**2/(4*(b*r**2/2+t**2))
    A2v = gx**2 * a/2 + b*t**2/2 + r**2 - c*t*gx + ((b-2)*t*r - c*gx*r)**2/(4*(b*r**2/2+t**2))
    # OLD  A2u = (-r**2/(4*A)) * ((b-2)*t + 2*c*gp)**2 + a/2 * gp**2 + A + 2*c*t*gp
    # OLD A2v = (-r**2/(4*A)) * ((b-2)*t + 2*c*gx)**2 + a/2 * gx**2 + A + 2*c*t*gx

    A3u = r**2/A + (t**2/(4*A**3)) * ((b-2)*r*t - c*gp*r)**2 + (r*t*((b-2)*r*t - c*gp*r))/(A**2)
    A3v = r**2/A + (t**2/(4*A**3)) * ((b-2)*r*t - c*gx*r)**2 + (r*t*((b-2)*r*t - c*gx*r))/(A**2)
    # OLD A3u = r**2/A + (t**2/(4*A**3)) * ((b-2)*t*r + 2*c*gp*r) + (t*r/A**(7/6))*((b-2)*t*r + 2*c*gp*r)
    # OLD A3v = r**2/A + (t**2/(4*A**3)) * ((b-2)*t*r + 2*c*gx*r) + (t*r/A**(7/6))*((b-2)*t*r + 2*c*gx*r)

    A4u = (gp**2 + 1) / 2
    A4v = (gx**2 + 1) / 2

    A5u = A4u + A2u
    A5v = A4v + A2v

    B2u = -2 * np.imag(alpha) * (-gp + 1)
    B2v = 2 * np.real(alpha) * (-gx + 1)

    C = 1/A - (t/A)**2
    E = np.exp(-B2u**2/(4*A5u) - B2v**2/(4*A5v))

    #######
    # print('A:', A)
    # print('A2u:', A2u)
    # print('A2v:', A2v)
    # print('A3u:', A3u)
    # print('A3v:', A3v)
    # print('A5u:', A5u)
    # print('A5v:', A5v)
    # print('B2u:', B2u)
    # print('B2v:', B2v)
    # print('C:', C)
    # print('E:', E)
    ########

    res = (C/np.sqrt(A5u * A5v) - A3u * (2 * A5u - B2u**2)/np.sqrt(2 * A5u * A5v) - A3v * (2 * A5v - B2v**2)/np.sqrt(2 * A5v * A5u) ) * E
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
    print(res)
    if not res['success']:
        raise AssertionError('Failure in optimization')
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
