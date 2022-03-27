import matplotlib.pyplot as plt
from numpy import inf, exp, sqrt, linspace, array, pi
import numpy as np
from scipy.integrate import quad
from decimal import Decimal

A = 45 / (4 * pi ** 4)
g_x = 4
g_s = 106.75

x_min = 0.1
x_max = 50

x = linspace(x_min, x_max, 200)


def N_eq(x):
    if isinstance(x, np.float64) or type(x) == float:
        return A * (g_x / g_s) * quad(lambda a: a ** 2 / (exp(sqrt(a ** 2 + x ** 2)) - 1), 0, inf, epsabs=inf)[0]
    else:
        return array(
            [A * (g_x / g_s) * quad(lambda a: a ** 2 / (exp(sqrt(a ** 2 + x_ ** 2)) - 1), 0, inf, epsabs=inf)[0] for x_
             in
             x])


def dN_eq(x):
    # x = logx
    if isinstance(x, np.float64) or type(x) == float:
        return A * (g_x / g_s) * quad(lambda a: -a ** 2 * x * exp(sqrt(a ** 2 + x ** 2)) / (
                sqrt(a ** 2 + x ** 2) * (exp(sqrt(a ** 2 + x ** 2)) - 1) ** 2), 0, inf, epsabs=inf)[0]
    else:
        return array(
            [A * (g_x / g_s) * quad(lambda a: -a ** 2 * x_ * exp(sqrt(a ** 2 + x_ ** 2)) / (
                    sqrt(a ** 2 + x_ ** 2) * (exp(sqrt(a ** 2 + x_ ** 2)) - 1) ** 2), 0, inf, epsabs=inf)[0] for x_ in
             x])


y = N_eq(x)
# y_deriv = dN_eq(x)
plt.loglog(x, y, 'k--', linewidth=0.5, label=r'$N_{eq}(x)$')
# plt.loglog(x, y_deriv, 'r--', linewidth=0.5, label=r'$dN_{eq}(x)$')
#
# plt.title(r'$N_x^{eq}(x)$')
# plt.ylabel(r'$N_x^{eq}$')
# plt.xlabel(r'$x$')
# plt.xlim(x_min, x_max)
# plt.ylim(1e-10, 1e-1)
# plt.savefig('N_eq_plot', dpi=300)
# plt.show()

x_vals = [0.1, 1., 10.]
for x in x_vals:
    print("N_eq({}) = {}".format(x, "{:.3e}".format(N_eq(x))))
