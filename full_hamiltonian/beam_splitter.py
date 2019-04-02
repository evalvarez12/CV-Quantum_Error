# -*- coding: utf-8 -*-
"""
Functions to implement beam splitter operations

Created on Mon Mar  4 11:58:35 2019

@author: Eduardo Villasenor
"""

import qutip as qt
import numpy as np
import tools


def beam_splitter(operators, theta, tensor=True):
    """
    Returns the corresponding beam splitter operator.
    Asummes symmetric beam splitter

    Returns
    --------
    bs_ops : qobj - Operator corresponding to the beam splitter transform
    """
#    U = np.array([[np.cos(theta) + 1j*np.sin(theta)] , [1j*np.sin(theta) + np.cos(theta)]])

    a, b = operators

    if tensor:
        # Tensoring additional identities required to form the new Hlbert space
        N = a.dims[0][0]
        identity = qt.identity(N)
        a = qt.tensor(a, identity)
        b = qt.tensor(identity, b)

    operators_out = [a*np.cos(theta) + b*1j*np.sin(theta), a*1j*np.sin(theta) + b*np.cos(theta)]
    return operators_out


def tritter(operators, theta1, theta2, tensor=True):
    """
           |
    c ->---/----->
           |
           |
    b ->---/----->
           |
           ^
           a
    """
    a, b, c = operators
    a, b = beam_splitter([a, b], theta1, tensor)

    if tensor:
        # Tensor additional identityes reqiured
        c = qt.tensor(qt.identity(c.dims[0][0]), c)

    a, c = beam_splitter([a, c], theta2, tensor)

    if tensor:
        b = qt.tensor(b, qt.identity(b.dims[0][0]))

    return [a, b, c]


def beam_splitter_U(N, theta, pos=[0,1], N_modes=2):
    """
             |
    b_op ->--/----->
             |
             ^
            a_op
    """
    a_op = qt.tensor(qt.destroy(N), qt.identity(N))
    b_op = qt.tensor(qt.identity(N), qt.destroy(N))
#    a_op = tools.tensor(qt.destroy(N), N, pos[0], N_modes)
#    b_op = tools.tensor(qt.destroy(N), N, pos[1], N_modes)


    # Define the mode-mixing Hamiltonian
    H = (a_op.dag()*b_op + a_op*b_op.dag())

    # Unitary evolution
    U = (1j*theta*H).expm()

    # If any append extra required H-spaces and permute as required
    if N_modes > 2:
        U = qt.tensor([U] + [qt.qeye(N)]*(N_modes-2))


        # TODO improve this list to permute
        permute_list = list(range(N_modes))

        permute_list[pos[0]] = 0
        permute_list[0] = pos[0]

        permute_list[pos[1]] = permute_list[1]
        permute_list[1] = pos[1]
        U = U.permute(permute_list)
    return U


# def beam_splitter_applyU(rho, theta, pos=[0,1], N_modes=2):
#
#     # Get unitary opertator
#     U = beam_splitter_Uoperator(rho.dims[0][0], theta, pos, N_modes)
#
#     rho = U * rho * U.dag()
#     return rho


def tritter_applyU(rho, theta1, theta2, pos=[0, 1, 2]):
    """
           |
    c ->---/----->
           |
           |
    b ->---/----->
           |
           ^
           a
    """
    N = rho.dims[0][0]
    N_modes = len(rho.dims[0])

    pos_a, pos_b, pos_c = pos

    U1 = beam_splitter_Uoperator(N, theta1, pos=[pos_a, pos_b], N_modes=N_modes)
#    U1 = qt.tensor(U1, qt.identity(N))
    U2 = beam_splitter_Uoperator(N, theta2, pos=[pos_a, pos_c], N_modes=N_modes)
#    U2 = qt.tensor(qt.identity(N), U2)

    rho = U1 * rho * U1.dag()
    rho = U2 * rho * U2.dag()
    return rho


def loss_channel(eta, a_in):
    theta = np.arccos(np.sqrt(eta))

    a_noise =qt.destroy(a_in.dims[0][0])

    a_out, _ = beam_splitter([a_in, a_noise], theta)


def loss_channel_applyU(rho_in, pos, eta):
    N = rho_in.dims[0][0]
    N_modes = len(rho_in.dims[0])

    theta = np.arccos(np.sqrt(eta))

    vacuum = qt.basis(N) * qt.basis(N).dag()

    rho = qt.tensor(rho_in, vacuum)

    rho = beam_splitter_applyU(rho, theta, [pos, N_modes], N_modes+1)

    rho = rho.ptrace(range(N_modes))

    return rho
