import warnings

from q1_b import *
# https://www.diva-portal.org/smash/get/diva2:720978/FULLTEXT01.pdf

import matplotlib.pyplot as plt
from numpy import sqrt, arange, pi, log
from scipy.integrate import odeint, ode


def l(g_w):
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

        return ((pi / 3) * sqrt(g_s / 10) * m_x ** 2) / (M_pl)

    # return ((2 * pi ** 2) / 45) * g_s * (m_x ** 3 / H()) * sigma_v()
    # https://www.ippp.dur.ac.uk/~dcerdeno/Dark_Matter_Lab_files/DM.pdf (2.39)
    return ((2 * pi ** 2) / 45) * g_s * m_x ** 3 * sigma_v() / H()


def dfX(logx, logN, g_w):
    return l(g_w) * (exp(logN) - exp(-logN + 2 * log(N_eq(exp(logx))))) * exp(-logx)


def dfN(logx, logN, g_w):
    return - l(g_w) * (exp(logN) + exp(-logN + 2 * log(N_eq(exp(logx))))) * exp(-logx)


def f(logx, logN):
    return -l(g_w) * (exp(logN) - exp(-logN + 2 * log(N_eq(exp(logx))))) * exp(-logx)


def dN(logx, logN, dlogx, g_w):
    return (f(logx, logN) * dlogx + dfX(logx, logN, g_w) * (dlogx ** 2)) / (1 - dfN(logx, logN, g_w) * dlogx)


g_w_array = [1e-4, 1e-3, 1e-2, 1e-1, 1]
dx = 0.05
x_min = 0.1
x_max = 1000.
x_evals = log(np.arange(x_min, x_max + dx, step=dx))
for g_w in g_w_array:
    log_N_x = []
    val = log(N_eq(exp(x_min)))
    for i in range(len(x_evals)):
        if i == 0:
            log_N_x.append([val])
        else:
            old_val = val
            val += dN(x_evals[i], old_val, x_evals[i] - x_evals[i - 1], g_w)
            log_N_x.append([val])
    if g_w >= 1e-3:
        print("Experimental N(oo) = {} @ gw = {}".format(exp(log_N_x[-1]), g_w))
        plt.hlines(6e-17 * g_w ** (-3.8), x_min, x_max, linestyles='dashed', colors='k', linewidth=0.5)
        print("Analyic N(oo) = {}".format(6e-17 * g_w ** (-3.8)))
    plt.loglog(exp(x_evals), exp(log_N_x), linewidth=1.5, label=r'$g_w$ = {}'.format(g_w))

plt.title(r'$N_x(x)$')
plt.ylabel(r'$N_x$')
plt.xlabel(r'$x$')
plt.xlim(x_min, x_max)
plt.legend(loc='best')
plt.savefig('g_w_transitions', dpi=300)
plt.show()
