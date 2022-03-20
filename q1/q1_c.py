from q1_b import *

import matplotlib.pyplot as plt
from numpy import sqrt, arange, pi, log
from scipy.integrate import odeint, solve_ivp

x_min = 0.1
x_max = 1000.
dx = 0.01
x = arange(x_min, x_max, dx)

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


for g_w_ in g_w:
    soln = odeint(f, y0=W_eq(x_min), t=x, tfirst=True)
    # soln2 = solve_ivp(F, [x_min, x_max], [W_eq(x_min)], t_eval=x)
    plt.loglog(soln, '-')
    # plt.xlim(min(x), max(x))
    plt.title("{}".format(g_w_))
plt.show()
