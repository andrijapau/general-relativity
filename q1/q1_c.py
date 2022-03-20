from q1_b import *

import matplotlib.pyplot as plt
from numpy import sqrt, arange, pi, log
from scipy.integrate import odeint

g_w = [1e-4, 1e-3, 1e-2, 1e-1, 1]


def W_eq(x):
    return log(N_eq(x))


def f(x, W):
    m_x = 500
    M_pl = 2.4e18
    g_s = 106.75
    H = (pi / 3) * sqrt(g_s / 10) * m_x / (M_pl * x)
    m_w = 80.4
    sigma_v = ((g_w_ * m_x) / (4 * pi * m_w)) ** 2
    A = (2 * pi ** 2 / 45) * g_s * (m_x ** 3 * sigma_v / H)
    f_result = (A / x ** 2) * (exp(2 * W_eq(x) - W) - W)
    return f_result


x_min = 0.1
x_max = 1000.
dx = 0.01
x_eval = arange(x_min, x_max, dx)

for g_w_ in g_w:
    soln = odeint(f, y0=W_eq(x_min), t=x_eval, tfirst=True)
    plt.loglog(x_eval, soln, '-')
    plt.xlim(min(x_eval), max(x_eval))
    plt.title("{}".format(g_w_))
plt.show()
