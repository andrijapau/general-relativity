import matplotlib.pyplot as plt
from numpy import inf, exp, sqrt, linspace, array, pi
import numpy as np
from scipy.integrate import quad

A = 45 / (4 * pi ** 4)
g_x = 4
g_s = 106.75

x_min = 0.1
x_max = 50

x = linspace(x_min, x_max, 200)


def N_eq(x):
    if isinstance(x, np.float64) or type(x) == float:
        return A * (g_x / g_s) * quad(lambda a: a ** 2 / sqrt(exp(a ** 2 + x ** 2) + 1), 0, inf, epsabs=inf)[0]
    else:
        return array(
            [A * (g_x / g_s) * quad(lambda a: a ** 2 / sqrt(exp(a ** 2 + x_ ** 2) + 1), 0, inf, epsabs=inf)[0] for x_ in
             x])


y = N_eq(x)
plt.loglog(x, y, '-')
plt.xlim(x_min, x_max)
plt.ylim(1e-10, 1e-1)
plt.show()

x_vals = [0.1, 1., 10.]
for x in x_vals:
    print("N_eq({}) = {}".format(x, N_eq(x)))
