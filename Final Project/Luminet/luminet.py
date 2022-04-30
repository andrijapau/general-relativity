import numpy as np
import matplotlib as plt
from luminet_equations import *
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

M = 1

r, phi = 10 * M, np.arange(0, 2 * np.pi, 0.001)

alpha, theta0 = 0 * np.pi / 180, 80 * np.pi / 180


def func(P):
    return equation13(P, M, alpha, theta0, r)


P_guess = np.arange(0, 40 * M)
plt.plot(P_guess, func(P_guess))
plt.show()
