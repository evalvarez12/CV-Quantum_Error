import numpy as np

def fidelity(T, eps, alpha):
    E = np.exp(-2 * (1 - np.sqrt(T))**2 * np.abs(alpha)**2 / (eps + 2))
    return  2 * E / (eps + 2)
