import numpy as np
import scipy as sp
from scipy.special import ellipj, ellipkinc


def u1(P, M):
    return -(Q(P, M) - P + 2 * M) / (4 * M * P)


def u2(P):
    return 1 / P


def u3(P, M):
    return (Q(P, M) + P - 2 * M) / (4 * M * P)


def Q(P, M):
    return np.sqrt((P - 2 * M) * (P + 6 * M))


def b(P, M):
    return P ** 3 / (P - 2 * M)


def k(P, M):
    return np.sqrt((u2(P) - u1(P, M)) / (u3(P, M) - u1(P, M)))
    # return np.sqrt((Q - P + 6 * M) / (2 * Q))


def zeta_inf(P, M):
    arg = (Q(P, M) - P + 2 * M) / (Q(P, M) - P + 6 * M)
    return np.arcsin(np.sqrt(arg))


def zeta_r(Q, P, M, r):
    return np.arcsin(np.sqrt((Q - P + 2 * M + (4 * M * P) / r) / (Q - P + 6 * M)))


def gamma(alpha, theta0):
    num = np.cos(alpha)
    denom = np.sqrt(np.cos(alpha) ** 2 + (1 / np.tan(theta0)) ** 2)
    return np.arccos(num / denom)


def F(zeta, m):
    return ellipkinc(zeta, m)


def equation13(P, M, alpha, theta0, r):
    m_ = k(P, M)
    sn_arg = (gamma(alpha, theta0) / 2) * np.sqrt(2 * M * (u3(P, M) - u1(P, M))) + F(zeta_inf(P, M), m_)
    RHS = u1(P, M) + (u2(P) - u1(P, M)) * (ellipj(sn_arg, m_)[0]) ** 2
    eqn = 1 - r * RHS
    return eqn


def r_eqn(Q, P, M, phi):
    alpha = np.sqrt((u3(Q, P, M) - u1(Q, P, M)) / 2)
    beta = (u2(P) - u3(Q, P, M)) / (u3(Q, P, M) - u1(Q, P, M))
    return 1 / (u1(Q, P, M) + (u2(P) - u1(Q, P, M)) * ellipj(alpha * phi, beta)[0] ** 2)
