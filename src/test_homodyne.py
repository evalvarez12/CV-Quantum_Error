import cv_system as cv
import wigner_plots as wp
N = 20
a = 1j

sys = cv.System(N, 1)
sys.apply_SMD(a)

# wp.plot(sys.state, [-5, 5])


# sys.apply_photon_addition(t, 1, 1)

state, prob = sys.homodyne_measurement()

# wp.plot(state, [-5, 5])
