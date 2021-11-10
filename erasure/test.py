import src.cv_system as cv


N = 20
r = 2
t = 0.1

sys = cv.System(N, 2)
sys.apply_TMS(r)

# sys.apply_photon_addition(t, 1, 1)

print(sys.get_full_CM())
