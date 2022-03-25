import warnings

from q1_b import *
# https://www.diva-portal.org/smash/get/diva2:720978/FULLTEXT01.pdf

import matplotlib.pyplot as plt
from numpy import sqrt, arange, pi, log
from scipy.integrate import odeint, ode

g_w = 1e-3


def l():
    m_x = 500
    g_s = 106.75

    M_pl = 2.4e18

    def sigma_v():
        m_x = 500
        m_w = 80.4

        return ((g_w * m_x) / (4 * pi * m_w)) ** 2

    def H():
        M_pl = 2.4e18
        m_x = 500
        g_s = 106.75

        return ((pi / 3) * sqrt(g_s / 10) * m_x) / (M_pl * x)

    # return (2 * pi ** 2 / 45) * g_s * (m_x ** 3 * sigma_v() / H())
    # https://www.ippp.dur.ac.uk/~dcerdeno/Dark_Matter_Lab_files/DM.pdf (2.39)
    return 0.26 * sqrt(g_s) * M_pl * m_x * sigma_v()


def dfX(logx, logN):
    return l() * (exp(logN) - exp(-logN + 2 * log(N_eq(exp(logx))))) * exp(-logx)


def dfN(logx, logN):
    return - l() * (exp(logN) + exp(-logN + 2 * log(N_eq(exp(logx))))) * exp(-logx)


def f(logx, logN):
    return -l() * (exp(logN) - exp(-logN + 2 * log(N_eq(exp(logx))))) * exp(-logx)


def dN(logx, logN, dlogx):
    return (f(logx, logN) * dlogx + dfN(logx, logN) * (dlogx ** 2)) / (1 - dfN(logx, logN) * dlogx)


dx = 0.05
x_min = 0.1
x_max = 1000.
x_evals = log(np.arange(x_min, x_max, step=dx))
log_N_x = []
val = log(N_eq(x_min))
for i in range(len(x_evals)):
    if i == 0:
        pass
    old_val = val
    val += dN(x_evals[i], old_val, x_evals[i] - x_evals[i - 1])
    log_N_x.append([val])
plt.loglog(exp(x_evals), exp(log_N_x))

plt.show()
