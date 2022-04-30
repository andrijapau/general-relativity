import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib
import scipy.integrate as integrate

plt.rcParams.update({'font.size': 16})

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


def G(M, b, u):
    return (2 * M * u ** 3 - u ** 2 + 1 / b ** 2)


def dudphi(phi, u, M, b):
    return np.sqrt(2 * M * G(M, b, u))


def dudphi_new(phi, u, M, b):
    return 1 / np.sqrt(2 * M * G(M, b, u))


M = 1
bc = np.sqrt(27) * M
us = np.linspace(0, 1, 1000)

bs = [bc - 1, bc, bc + 2]
plt.plot(us, G(M, bs[0], us), 'k:', linewidth=1.5, label=r'$b < b_c$')
plt.plot(us, G(M, bs[1], us), 'k-', linewidth=1.5, label=r'$b = b_c$')
plt.plot(us, G(M, bs[2], us), 'k--', linewidth=1.5, label=r'$b > b_c$')
plt.locator_params(axis='y', nbins=5)
plt.locator_params(axis='x', nbins=5)
plt.gca().set_aspect('equal')
plt.ylim(-0.1, 0.1)
plt.xlim(0, 0.2)
plt.title(r"G(u) Potential", fontsize=16)
plt.ylabel(r'$G(u)$', rotation=0, fontsize=16, labelpad=10)
plt.xlabel(r'$u$', fontsize=16)
plt.axhline(0, color='black', linewidth=0.75)
plt.axvline(0, color='black', linewidth=0.75)
plt.legend(loc='best')
plt.savefig('g-potential.png', bbox_inches='tight', dpi=400)
plt.show()
plt.show()


# P = 3.2 * M
# b = np.sqrt(P ** 3 / (P - 2 * M))
# bc = np.sqrt(27) * M
# print(b, b > bc)
# Q = np.sqrt((P - 2 * M) * (P + 6 * M))
# u3 = (Q + P - 2 * M) / (4 * M * P)
# phi_0, u_0 = 0, 0
# u_values = np.arange(u_0, 1 / P, 0.00001)
# print(u3)
# u_values_2 = np.arange(1 / P, u3, 0.001)


def rungekutta4(f, y0, t, args=()):
    n = len(t)
    y = np.zeros((n, len(y0)))
    y[0] = y0
    for i in range(n - 1):
        h = t[i + 1] - t[i]
        k1 = f(y[i], t[i], *args)
        k2 = f(y[i] + k1 * h / 2., t[i] + h / 2., *args)
        k3 = f(y[i] + k2 * h / 2., t[i] + h / 2., *args)
        k4 = f(y[i] + k3 * h, t[i] + h, *args)
        y[i + 1] = y[i] + (h / 6.) * (k1 + 2 * k2 + 2 * k3 + k4)
    return y


P = 3.2 * M
P_vals = [val * M for val in range(3, 13, 2)]
# P_vals = [val * M for val in range(3, 13, 2)]
for P in P_vals:
    phi = []
    phi_new = []
    phi_total = []

    b = np.sqrt(P ** 3 / (P - 2 * M))
    # b = np.sqrt(30) * M
    phi_0, u_0 = 0, 0
    u_values = np.arange(u_0, 1 / P, 0.00001)
    soln = rungekutta4(dudphi_new, [phi_0], u_values, args=(M, b))
    phi = soln + phi_0
    # plt.plot(phi, 1 / u_values)
    # plt.show()
    max, min = phi[-1], phi[0]
    a = (max + min) / 2
    phi_new = [a - (phi_ - a) for phi_ in phi]
    phi_new += phi[-1][0]
    for phi_ in phi:
        phi_total.append(phi_[0])
    for phi_ in phi_new:
        phi_total.append(phi_[0])
    u_total = np.tile(u_values, 2)
    r = 1 / u_total
    y = r * np.sin(phi_total)
    x = r * np.cos(phi_total)
    if P > 3 * M:
        y = 2 * b - r * np.sin(phi_total)
    plt.plot(x, y, 'k-', linewidth=1.5)

photon_sphere = plt.Circle((0, 0), 3, facecolor='none', edgecolor='black', linestyle='dashed')
bh = plt.Circle((0, 0), 2, facecolor='black', edgecolor='black')
plt.gca().add_patch(photon_sphere)
plt.gca().add_patch(bh)
scale = 25
plt.ylim(-scale, scale)
plt.xlim(-scale, scale)
# plt.axline((0, 0), xy2=(10, 0), color='black', linewidth=0.75)
# plt.axline((0, 0), slope=0.85, color='k', linewidth=0.75)
plt.locator_params(axis='y', nbins=10)
plt.locator_params(axis='x', nbins=10)
plt.gca().set_aspect('equal')
plt.title(r"Photon Orbit", fontsize=16)
plt.ylabel(r'$y/M$', rotation=0, fontsize=16, labelpad=5)
plt.xlabel(r'$x/M$', fontsize=16)
plt.savefig('photon-orbit-multiple.png', bbox_inches='tight', dpi=400)
plt.show()
