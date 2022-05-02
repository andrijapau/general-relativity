# Import Modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.integrate import odeint

plt.rcParams.update({'font.size': 16})
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Define global variables
global rs, M, rc
M = 1
rs = 2 * M
rc = (3 / 2) * rs


def Q(P, M):
    return np.sqrt((P - 2 * M) * (P + 6 * M))


def u1_root(P):
    return -(Q(P, M) - P + 2 * M) / (4 * M * P)


def u2_root(P):
    return 1 / P


def u3_root(P):
    return (Q(P, M) + P - 2 * M) / (4 * M * P)


def b(P):
    return np.sqrt(P ** 3 / (P - 2 * M))


def circular_orbit(phi, r0_photon):
    u_photon = 1 / r0_photon

    uc = 1 / rc
    u1 = 1 / rs - 2 * uc

    a0 = - np.sqrt(2 / (M * (uc - u1)))
    phi0 = a0 * np.arctanh(np.sqrt((u1 - u_photon) / (u1 - uc)))

    return 1 / (u1 + (uc - u1) * np.tanh((phi - phi0) / a0) ** 2)


def psi(u, phi, P, u1, u2, u3):
    return 1 / np.sqrt(rs * (u - u1) * (u - u2) * (u - u3))


def unbounded_orbit(P, u_values, phi0, u1, u2, u3):
    soln = odeint(psi, y0=phi0, t=u_values, tfirst=True, args=(P, u1, u2, u3))
    return soln[:, 0]


# Plot circular orbit
P = 7 * rs
u1, u2, u3 = u1_root(P), u2_root(P), u3_root(P)
r0_photon = 15 * rs
phi = np.arange(0, 4 * np.pi, 0.01)
r = circular_orbit(phi, r0_photon)
x = r * np.cos(phi)
y = r * np.sin(phi)
photon_sphere = plt.Circle((0, 0), 3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=3)
bh = plt.Circle((0, 0), 2, facecolor='black', edgecolor='black')
plt.gca().add_patch(photon_sphere)
plt.gca().add_patch(bh)
scale = 5
plt.ylim(-scale, scale)
plt.xlim(-scale, scale)
plt.locator_params(axis='y', nbins=10)
plt.locator_params(axis='x', nbins=10)
plt.gca().set_aspect('equal')
plt.title(r"Photon Orbit", fontsize=16)
plt.ylabel(r'$y/M$', rotation=0, fontsize=16, labelpad=5)
plt.xlabel(r'$x/M$', fontsize=16)
plt.plot(x, y, 'k-', linewidth=0.75)
plt.savefig('Photos/photon-orbit-circle.png', bbox_inches='tight', dpi=400)
plt.show()

# Plot multiple orbits
num_of_orbits = 10
dec_values = np.linspace(3.05, 20.25, num_of_orbits)
P_values = [val * M for val in dec_values]
for P in P_values:
    phi_total = []

    u1, u2, u3 = u1_root(P), u2_root(P), u3_root(P)
    r0_photon = 15 * rs

    u_values = np.arange(1 / r0_photon, u2, 0.00001)
    phi0 = 0
    phi = unbounded_orbit(P, u_values, phi0, u1, u2, u3) + phi0

    r = 1 / u_values

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    # color = next(plt.gca()._get_lines.prop_cycler)['color']
    color = 'black'

    plt.plot(x, y, color=color, linewidth=1.5)

    rot_angle = phi[-1] - np.pi
    a = np.sin(rot_angle) / np.cos(rot_angle)
    b = 0
    # x_values = np.arange(r[-1] * np.cos(phi[-1]), 50, 0.1)
    # plt.plot(x_values, a * x_values + b, color=color, linestyle='dashed', linewidth=0.75)

    x_ = x * (1 - a ** 2) / (1 + a ** 2) + (y - b) * (2 * a) / (a ** 2 + 1)
    y_ = x * (2 * a) / (a ** 2 + 1) + (y - b) * (a ** 2 - 1) / (a ** 2 + 1) + b
    plt.plot(x_, y_, color=color, linestyle='-', linewidth=1.5, label=r'$P = {}M$'.format(round(P, 2)))

photon_sphere = plt.Circle((0, 0), 3, facecolor='none', edgecolor='black', linestyle='dashed', linewidth=2)
bh = plt.Circle((0, 0), 2, facecolor='black', edgecolor='black')
plt.gca().add_patch(photon_sphere)
plt.gca().add_patch(bh)
scale = 15
plt.ylim(-scale, scale)
plt.xlim(-scale, scale)
plt.locator_params(axis='y', nbins=10)
plt.locator_params(axis='x', nbins=10)
plt.gca().set_aspect('equal')
plt.title(r"Photon Orbit", fontsize=16)
plt.ylabel(r'$y/M$', rotation=0, fontsize=16, labelpad=5)
plt.xlabel(r'$x/M$', fontsize=16)
# plt.legend(loc='best', fontsize=10)
plt.savefig('Photos/photon-orbit-multiple.png', bbox_inches='tight', dpi=400)
plt.show()
