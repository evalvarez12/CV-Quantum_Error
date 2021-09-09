"""
Functions to perform homodyne measurements
"""

import qutip as qt
import numpy as np
import tools


def get_homodyne_projector(N, q_m, theta, eta):
    a = qt.destroy(N)

    # Using SNU units
    q = a.dag() + a
    p = 1j*(a.dag() - a)
    qtheta = q*np.cos(theta) + p*np.sin(theta)
    # Projector obtained from : 10.1103/RevModPhys.81.299
    P = (eta * (-(q_m/eta)**2 + 2*q_m*qtheta/eta - qtheta**2)).expm()
    return P/np.sqrt(np.pi/eta)


def homodyne_measurement(state, mode=0, theta=0, eta=1, mlim=[-3, 3], sample_size=20):
    """
    Perfomr a homodyne measurement
    """
    N = state.dims[0][0]
    Nmodes = len(state.dims[0])

    q_samples = np.linspace(mlim[0], mlim[1], sample_size)
    samples = np.zeros_like(q_samples)

    for i in range(q_samples):
        projector = get_homodyne_projector(N, q_samples[i], theta, eta)
        projector = tools.tensor(N, projector, mode, Nmodes)
        samples[i] = qt.expect(projector, state)

    # Perform 1D polyfit to approximate the PDF
    fit = np.polyfit(q_samples, samples, 1)
    pdf_func = np.poly1d(fit)
    # New points to find pdf
    qs = np.linspace(mlim[0], mlim[1], 10 * sample_size)
    pdf = pdf_func(qs)

    # Select a measurement result
    q_measured = np.random.choice(qs, p=pdf)
    projector = get_homodyne_projector(N, q_measured, theta, eta)
    projector = tools.tensor(N, projector, mode, Nmodes)

    if state.isket:
        state = state * state.dag()

    state = projector.dag() * state * projector
    prob = (state * projector).tr()
    state = state / prob
    return state, prob


def homodyne_measurement_gaussian_approx(state, mode=0, theta=0, eta=1):
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
    projector = get_homodyne_projector(N, q_measured, theta, eta)
    projector = tools.tensor(N, projector, mode, Nmodes)

    if state.isket:
        state = state * state.dag()

    state = projector.dag() * state * projector
    prob = (state * projector).tr()
    state = state / prob
    return state, prob
