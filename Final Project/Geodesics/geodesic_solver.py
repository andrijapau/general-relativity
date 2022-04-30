import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from cmath import sqrt
from scipy.special import ellipj, ellipkinc
from scipy import ndimage

global r_s, M
M = 1
r_s = 2 * M


def zeta(r):
    return r_s / r


def a(r_obs, xi):
    return (zeta(r_obs) ** 2) * ((1 - zeta(r_obs)) / (np.sin(xi) ** 2))


def D(a):
    return (a ** 2) * (((a / 2) ** 2) - 1 / 27)


def q(a):
    return (a ** 2) / 2 - 1 / 27


def rho(a):
    return np.sign(q(a)) / 3


def psi(a):
    if D(a) <= 0:
        return np.arccos(q(a) / rho(a) ** 3)
    if D(a) > 0:
        return np.arccosh(q(a) / rho(a) ** 3)


def zeta_roots(a):
    if D(a) <= 0:
        zeta1 = rho(a) * (np.cos(psi(a) / 3) + sqrt(3) * np.sin(psi(a) / 3)) + 1 / 3
        zeta2 = rho(a) * (np.cos(psi(a) / 3) - sqrt(3) * np.sin(psi(a) / 3)) + 1 / 3
        zeta3 = - 2 * rho(a) * np.cos(psi(a) / 3) + 1 / 3
        return zeta1, zeta2, zeta3
    if D(a) > 0:
        zeta1 = rho(a) * (np.cosh(psi(a) / 3) + 1j * sqrt(3) * np.sinh(psi(a) / 3)) + 1 / 3
        zeta2 = rho(a) * (np.cosh(psi(a) / 3) - 1j * sqrt(3) * np.sinh(psi(a) / 3)) + 1 / 3
        zeta3 = - 2 * rho(a) * np.cosh(psi(a) / 3) + 1 / 3
        return zeta1, zeta2, zeta3


def zeta_eqn(a, phi, phi0):
    zeta1, zeta2, zeta3 = zeta_roots(a)
    # print(zeta1, zeta2, zeta3)

    m = np.sqrt(1 / np.real(-(zeta1 - zeta3) / -(zeta2 - zeta3)))
    # print("m2", m2)

    sn_arg = np.real(0.5 * np.sqrt(-(zeta2 - zeta3)) * (phi + phi0))
    # print("sn", sn_arg)

    denom = 1 - (ellipj(sn_arg, m)[0]) ** 2
    num = zeta2 - zeta1 * (ellipj(sn_arg, m)[0] ** 2)

    zeta_ = num / denom

    return zeta_


def phi(i, k):
    kd, kr, ku = k
    num = -ku * np.tan(i)
    denom = np.sqrt((kr ** 2) + (1 + np.tan(i) ** 2) * (ku ** 2))
    if denom == 0:
        print(True)
        return 1
    return np.arccos(num / denom)


def k(ph, pv, resh, resv, fov):
    # fov *= np.pi / 180
    rho = resh / resv
    return [1, rho * (2 * ph / resh - 1) * np.tan(fov / 2), (1 - 2 * pv / resv) * np.tan(fov / 2)]


def xi(k):
    kd, kr, ku = k

    return np.arccos(1 / sqrt(1 + kr ** 2 + ku ** 2))


def w(k):
    kd, kr, ku = k
    return np.arctan(ku / kr)


r_obs = 120 * r_s
resh = 600
resv = 500
fov = 6
image = np.zeros((resh, resv))
for h in range(resh):
    for v in range(resv):
        # print(h, v)
        k_vec = k(h, v, resh, resv, fov)
        i = -85 * np.pi / 180
        if h == 300 and v == 249:
            break

        rq = np.real(r_s / zeta_eqn(a=a(r_obs, xi(k_vec)), phi=phi(i, k_vec), phi0=0))

        Rin, Rout = 3 * r_s, 15 * r_s
        if rq < 3 * r_s:
            # photon sphere
            image[h, v] = 100
        if Rin < rq < Rout:
            image[h, v] = rq
            # print(yes)
        else:
            pass
# print(phiq)
# x = rq * np.cos(phiq)
# y = rq * np.sin(phiq)
# plt.scatter(x, y, label=i * 180 / np.pi)
# plt.axhline(0, color='black', linewidth=0.75)
# plt.axvline(0, color='black', linewidth=0.75)
# plt.legend(loc='best')
# plt.show()
goodimage = image[0:300, 0:500]
image = np.vstack((goodimage, np.flip(goodimage, axis=0)))
# plt.imshow(ndimage.rotate(image, 90), cmap=plt.cm.gray)
plt.imshow(ndimage.rotate(image, 90))

plt.show()
