import warnings

from q1_b import *
# https://www.diva-portal.org/smash/get/diva2:720978/FULLTEXT01.pdf

import matplotlib.pyplot as plt
from numpy import sqrt, arange, pi, log
from scipy.integrate import odeint, ode

g_w = [1]


def W_eq(x):
    return log(N_eq(x))


def sigma_v():
    m_x = 500
    m_w = 80.4

    return ((g_w_ * m_x) / (4 * pi * m_w)) ** 2


def H():
    M_pl = 2.4e18
    m_x = 500
    g_s = 106.75

    return ((pi / 3) * sqrt(g_s / 10) * m_x) / (M_pl * x)


def l():
    m_x = 500
    g_s = 106.75
    M_pl = 2.4e18

    # return (2 * pi ** 2 / 45) * g_s * (m_x ** 3 * sigma_v() / H())
    # https://www.ippp.dur.ac.uk/~dcerdeno/Dark_Matter_Lab_files/DM.pdf (2.39)
    return 0.26 * sqrt(g_s) * M_pl * m_x * sigma_v()


def f(x, N):
    # f_result = (l() / (x ** 2)) * (exp(2 * W_eq(x) - W) - W)
    f_result = -(l() / x ** 2) * (N ** 2 - N_eq(x) ** 2)
    return f_result


x_min = 0.1
x_max = 1000.
dx = 0.1
x_eval = arange(x_min, x_max, dx)

plt.loglog(x_eval, N_eq(x_eval))
for g_w_ in g_w:
    soln = odeint(f, y0=round(N_eq(x_min), 3), t=x_eval, tfirst=True)
    plt.loglog(x_eval, soln, '-')
    plt.xlim(x_min, max(x_eval))
    plt.ylim(1e-15, 1e-1)
    plt.title("{}".format(g_w_))
plt.show()
