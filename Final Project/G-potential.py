import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.rcParams.update({'font.size': 16})
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


def G(M, b, u):
    return (2 * M * u ** 3 - u ** 2 + 1 / b ** 2)


def dudphi(phi, u, M, b):
    return np.sqrt(2 * M * G(M, b, u))


def dudphi_new(phi, u, M, b):
    return 1 / np.sqrt(2 * M * G(M, b, u))


# Define constants
M = 1
bc = np.sqrt(27) * M
us = np.linspace(-1 / 2, 1, 1000)

# Plot!
bs = [bc - 1, bc, bc + 2]
plt.plot(us, G(M, bs[0], us), 'k:', linewidth=1.5, label=r'$b < b_c$')
plt.plot(us, G(M, bs[1], us), 'k-', linewidth=1.5, label=r'$b = b_c$')
plt.plot(us, G(M, bs[2], us), 'k--', linewidth=1.5, label=r'$b > b_c$')
plt.ylim(-0.1, 0.1)
plt.title(r"G(u) Potential", fontsize=16)
plt.ylabel(r'$G(u)$', rotation=0, fontsize=16, labelpad=10)
plt.xlabel(r'$u$', fontsize=16)
plt.axhline(0, color='black', linewidth=0.75)
plt.axvline(0, color='black', linewidth=0.75)
plt.legend(loc='best')
plt.savefig('Photos/g-potential.png', bbox_inches='tight', dpi=400)
plt.show()
plt.show()
