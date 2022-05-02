# Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.integrate import odeint

plt.rcParams.update({'font.size': 16})
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


def V(r, M):
    return (1 / r ** 2) * (1 - (2 * M) / r)


def drdphi(phi, r, M, b):
    if (1 / b ** 2) - V(r, M) > 0:
        return -r * np.sqrt(
            (1 / b ** 2) - V(r, M)
        )
    else:
        return 0


phi_min, phi_max = 0, 10 * np.pi
phi_vals = np.arange(phi_min, phi_max, 0.01)

M = 1
b_vals = [np.sqrt(val) * M for val in range(20, 70)]
r0 = 10 * M
for b in b_vals:
    soln = odeint(drdphi, y0=r0, t=phi_vals, tfirst=True, args=(M, b))
    r = soln[:, 0]
    plt.plot(phi_vals, r, linewidth=0.25, color='black')
plt.hlines(y=3 * M, xmin=phi_min, xmax=phi_max, colors='black', linewidth=2, linestyles={"dashed"},
           label=r'$r_\star = 3M$')
plt.locator_params(axis='y', nbins=5)
plt.locator_params(axis='x', nbins=5)
plt.xlim(0, 25)
plt.ylim(0, 10)
plt.title(r"Photon Trajectory", fontsize=16)
plt.ylabel(r'$y/M$', rotation=0, fontsize=16, labelpad=10)
plt.xlabel(r'$x/M$', fontsize=16)
plt.annotate("Event Horizon", xy=(0, 2.35 * M), fontsize=14)
plt.legend(loc='best')
plt.savefig('impact-parameters.png', bbox_inches='tight', dpi=400)
plt.show()
