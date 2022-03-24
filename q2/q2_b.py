from particle_info import *
import numpy as np
import matplotlib.pyplot as plt

quark = quark()
gluon = gluon()
boson = boson()
pion = pion()
lepton = lepton()
X = X()


def g_star(T):
    g_star_max = 106.75
    if T > X.T():
        return g_star_max

    if X.T() >= T > quark.top.T():
        return g_star_max - X.g()

    if quark.top.T() >= T > boson.W.T():
        return g_star_max - (X.g() + quark.top.g())

    if boson.W.T() >= T > quark.bottom.T():
        return g_star_max - (X.g() + quark.top.g() + boson.W.g())

    if quark.bottom.T() >= T > pion.pi_pm.T_join():
        return g_star_max - (X.g() + quark.top.g() + boson.W.g() + quark.bottom.g())

    if pion.pi_pm.T_join() >= T > pion.pi_pm.T_leave():
        g_minus = X.g() + quark.top.g() + boson.W.g() + quark.bottom.g() + quark.charm.g() + quark.strange.g() + quark.up.g() + quark.down.g() + gluon.g() + lepton.tau.g()
        g_plus = pion.pi_pm.g() + pion.pi_0.g()
        return g_star_max - g_minus + g_plus

    if pion.pi_pm.T_leave() >= T > lepton.e.T():
        g_minus = pion.pi_pm.g() + pion.pi_0.g() + X.g() + quark.top.g() + boson.W.g() + quark.bottom.g() + quark.charm.g() + quark.strange.g() + quark.up.g() + quark.down.g() + gluon.g() + lepton.tau.g()
        return g_star_max - g_minus

    if lepton.e.T() >= T:
        g_minus = lepton.e.g() + pion.pi_pm.g() + pion.pi_0.g() + X.g() + quark.top.g() + boson.W.g() + quark.bottom.g() + quark.charm.g() + quark.strange.g() + quark.up.g() + quark.down.g() + gluon.g() + lepton.tau.g()
        return g_star_max - g_minus

    else:
        return -1

# T_threshold = [100 * 1000, 30 * 1000, 15 * 1000, 1 * 1000, 200, 50, 25]  # MeV
#
# T = np.arange(25, 100 * 1000 + 1)[::-1]  # MeV
# g_star = [g_star(T_val) for T_val in T]
#
# for T_val in T_threshold:
#     plt.vlines(T_val, 0, 200, linestyles='dashed', colors='k')
# plt.loglog(T, g_star, 'k-')
# plt.show()
