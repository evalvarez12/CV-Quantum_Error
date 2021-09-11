import cv_system as cv
import wigner_plots as wp
import matplotlib.pyplot as plt
import homodyne as hom
import qutip as qt

N = 30
a = 3

sys = cv.System(N, 2)
sys.apply_TMS(a)

# wp.plot(sys.state, [-5, 5])


# sys.apply_photon_addition(t, 1, 1)

state, prob = sys.homodyne_measurement()


P1 = hom.pdf_projector(N, q_m=1, theta=0, eta=1)
P2 = hom.get_homodyne_projector_OLD(N, x_measured=1)
P3 = qt.displace(N, 1) * qt.squeeze(N, 15) * qt.basis(N)
P3 = P3 * P3.dag()

P1 = qt.tensor(P1, qt.qeye(N))
P3 = qt.tensor(P3, qt.qeye(N))

# wp.plot(P, [-5, 5])
# wp.plot(P2, [-5, 5])
# wp.plot(P3, [-5, 5])

s1 = P1 * sys.state
s3 = P3 * sys.state

wp.plot(sys.state.ptrace(0), [-5, 5])
wp.plot(s1.ptrace(0), [-5, 5])
wp.plot(s3.ptrace(0), [-5, 5])

# wp.plot(state, [-5, 5])
