"""
Functions to perform homodyne measurements
"""

import qutip as qt
import numpy as np
import tools

# OLD function save just in case


def get_homodyne_projector_OLD(N, x_measured):
    a = qt.create(N)
    x = (np.sqrt(2) * x_measured * a) - (x_measured**2 + a*a)/2
    x = 1/np.pi**(1/4) * x.expm()
    # x = (x * qt.basis(N)).dag()
    return x


def pdf_projector(N, q_m, theta, eta):
    a = qt.destroy(N)

    # Using SNU units
    q = (a.dag() + a)/np.sqrt(2)
    p = 1j*(a.dag() - a)/np.sqrt(2)
    qtheta = q*np.cos(theta) + p*np.sin(theta)
    # Projector obtained from : 10.1103/RevModPhys.81.299
    P = (eta * (-(q_m/eta)**2 + 2*q_m*qtheta/eta - qtheta**2)).expm()
    return P/np.sqrt(np.pi/eta)


def pdf_projector_red_dim(N, q_m, theta, eta):
    a = qt.destroy(N)

    # Using SNU units
    q = (a.dag() + a)/np.sqrt(2)
    p = 1j*(a.dag() - a)/np.sqrt(2)
    qtheta = q*np.cos(theta) + p*np.sin(theta)
    # Projector obtained from : 10.1103/RevModPhys.81.299
    P = (eta * (-(q_m/eta)**2 + 2*q_m*qtheta/eta - qtheta**2)).expm()
    P = P/np.sqrt(np.pi/eta)
    P = (P * qt.basis(N)).dag()
    return P


def homodyne_measurement(state, mode=0, theta=0, eta=1, mlim=[-5, 5], sample_size=20):
    N = state.dims[0][0]
    Nmodes = len(state.dims[0])

    q_samples = np.linspace(mlim[0], mlim[1], sample_size)
    samples = np.zeros_like(q_samples)
    for i in range(sample_size):
        projector = pdf_projector(N, q_samples[i], theta, eta)
        projector = tools.tensor(N, projector, mode, Nmodes)
        samples[i] = qt.expect(projector, state)

    # Check the area of the PDF
    print("SUM PDF:", np.sum(samples * (q_samples[1] - q_samples[0])))

    # Perform 1D polyfit to approximate the PDF
    fit = np.polyfit(q_samples, samples, 20)
    pdf_func = np.poly1d(fit)
    # New points to find pdf
    qs = np.linspace(mlim[0], mlim[1], 10 * sample_size)
    ds = np.abs(qs[0] - qs[1])
    pdf = pdf_func(qs)

    # import matplotlib.pyplot as plt
    # plt.plot(q_samples, samples)
    # plt.plot(qs, pdf)
    # plt.show()
    # # from IPython import embed
    # # embed()

    pdf = pdf * ds
    pdf[pdf < 0] = 0
    pdf = pdf / np.sum(pdf)

    # Select a measurement result
    q_measured = np.random.choice(qs, p=pdf)
    # q_measured = np.random.choice(
    #     q_samples, p=samples/sum(samples))
    projector = pdf_projector(N, q_measured, theta, eta)
    projector = tools.tensor(N, projector, mode, Nmodes)

    if state.isket:
        state = state * state.dag()

    prob = (state * projector).tr()
    state = projector.dag() * state * projector
    # state = state / prob
    return state, prob


def homodyne_measurement_gauss_approx(state, mode=0, theta=0, eta=1):
    N = state.dims[0][0]
    Nmodes = len(state.dims[0])

    a = qt.destroy(N)
    q = a.dag() + a
    p = 1j*(a.dag() - a)
    qtheta = q*np.cos(theta) + p*np.sin(theta)
    qtheta = tools.tensor(N, qtheta, mode, Nmodes)

    # Calculate <qtheta> and var(qtheta)
    avg_q = qt.expect(qtheta, state)
    var_q = 2 * (qt.expect(qtheta*qtheta, state) - qt.expect(qtheta, state)**2)

    # Select random measurment result
    q_measured = np.random.normal(avg_q, np.sqrt(var_q))
    projector = pdf_projector(N, q_measured, theta, eta)
    projector = tools.tensor(N, projector, mode, Nmodes)

    if state.isket:
        state = state * state.dag()

    state = projector.dag() * state * projector
    prob = (state * projector).tr()
    state = state / prob
    return state, prob
