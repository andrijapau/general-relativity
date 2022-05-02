# Import modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from cmath import sqrt
from scipy.special import ellipj, ellipkinc
from scipy import ndimage

plt.rcParams.update({'font.size': 16})
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Define global constants
global r_s, M
M = 1
r_s = 2 * M


def zeta(r):
    return r_s / r


def a(r_obs, xi):
    return np.sqrt((zeta(r_obs) ** 2) * ((1 - zeta(r_obs)) / (np.sin(xi) ** 2)))


def D(a):
    return (a ** 2) * (a ** 2 / 4 - 1 / 27)


def q(a):
    return (a ** 2) / 4 - 1 / 27


def rho(a):
    return np.sign(q(a)) / 3


def psi(a):
    if D(a) <= 0:
        return np.arccos(q(a) / (rho(a) ** 3))
    if D(a) > 0:
        return np.arccosh(q(a) / (rho(a) ** 3))


def zeta_roots(a):
    if D(a) <= 0:
        zeta1 = rho(a) * (np.cos(psi(a) / 3) + np.sqrt(3) * np.sin(psi(a) / 3)) + 1 / 3
        zeta2 = rho(a) * (np.cos(psi(a) / 3) - np.sqrt(3) * np.sin(psi(a) / 3)) + 1 / 3
        zeta3 = - 2 * rho(a) * np.cos(psi(a) / 3) + 1 / 3
        return zeta1, zeta2, zeta3
    if D(a) > 0:
        zeta1 = rho(a) * (np.cosh(psi(a) / 3) + 1j * sqrt(3) * np.sinh(psi(a) / 3)) + 1 / 3
        zeta2 = rho(a) * (np.cosh(psi(a) / 3) - 1j * sqrt(3) * np.sinh(psi(a) / 3)) + 1 / 3
        zeta3 = - 2 * rho(a) * np.cosh(psi(a) / 3) + 1 / 3
        return zeta1, zeta2, zeta3


def zeta_eqn(a, phi, phi0):
    zeta1, zeta2, zeta3 = zeta_roots(a)

    m = np.real(np.sqrt((zeta2 - zeta1) / (zeta3 - zeta1)))

    sn_arg = np.real(sqrt(0.5 * (zeta2 - zeta1)) * (phi + phi0))

    denom = (ellipj(sn_arg, m)[1]) ** 2
    num = zeta2 - zeta1 * (ellipj(sn_arg, m)[0]) ** 2

    zeta_ = num / denom

    return zeta_


def phi(i, k):
    kd, kr, ku = k
    num = -ku * np.tan(i)
    denom = np.sqrt((kr ** 2) + (1 + np.tan(i) ** 2) * (ku ** 2))
    return np.arccos(num / denom)


def k(ph, pv, resh, resv, fov):
    rho = resh / resv
    k = [1, rho * (2 * (ph / resh) - 1) * np.tan(fov / 2), (1 - 2 * (pv / resv)) * np.tan(fov / 2)]
    k_mag = np.sqrt(k[0] ** 2 + k[1] ** 2 + k[2] ** 2)
    return [(1 / k_mag), (1 / k_mag) * rho * (2 * (ph / resh) - 1) * np.tan(fov / 2),
            (1 / k_mag) * (1 - 2 * (pv / resv)) * np.tan(fov / 2)]


def xi(k):
    kd, kr, ku = k

    return np.arccos(1 / np.sqrt(1 + (kr ** 2) + (ku ** 2)))


def w(k):
    kd, kr, ku = k
    return np.arctan(ku / kr) + np.pi


def phi_zeta(zeta, a):
    zeta1, zeta2, zeta3 = zeta_roots(a)
    m = np.real(np.sqrt((zeta1 - zeta2) / (zeta2 - zeta3)))
    num = (zeta2 - zeta)
    denom = (zeta3 - zeta)
    sn_arg = np.real(np.sqrt(num / denom))
    return np.real(2 / (sqrt(-zeta3 + zeta1)) * ellipkinc(sn_arg, m))


factor = 2
r_obs = 120 * r_s
resh = 1080 // factor
resv = 1080 // factor
fov = 25 * np.pi / 180
i_deg = 85
image = np.zeros((resh, resv))
for h in range(resh):
    for v in range(resv):

        k_vec = k(h, v, resh, resv, fov)
        i = i_deg * np.pi / 180
        phiq = phi(i, k_vec)
        phi0 = 0
        rq = np.real(r_s / zeta_eqn(a=a(r_obs, xi(k_vec)), phi=phiq, phi0=phi0))
        Rs = [3 * r_s, 5 * r_s, 10 * r_s, 15 * r_s]

        if Rs[0] < rq <= Rs[1]:
            image[h, v] = rq * 1
        if Rs[1] < rq <= Rs[2]:
            image[h, v] = rq * 2
        if Rs[2] < rq <= Rs[3]:
            image[h, v] = rq * 3
        elif rq < Rs[0]:
            image[h, v] = 0

plt.imshow(ndimage.rotate(image, 90), cmap=plt.cm.gray)
plt.title(r"Schwarzschild Black Hole @ Inclination = {}".format(i_deg), fontsize=16)
plt.ylabel(r'$p_v$', rotation=0, fontsize=16, labelpad=5)
plt.xlabel(r'$p_h$', fontsize=16)
plt.savefig('Photos/my-black-hole-photo-{}.png'.format(i_deg), bbox_inches='tight', dpi=400)
plt.show()
