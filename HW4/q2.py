# Import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def system(current_state, t):
    '''
    System of equations that solves the coupled system of ODEs for delta and phi
    '''

    phi, delta, dphi, ddelta = current_state
    k = 1
    ddphi = -(3 / t) * dphi - k ** 2 * phi
    dddelta = 3 * ddphi - k ** 2 * phi - (1 / (2 * t)) * (ddelta - 3 * dphi)

    return [dphi, ddelta, ddphi, dddelta]


# Solve System
initial_state = [-3 / 4, -3 / 4, 0, 0]
ktau_vals = np.arange(1e-2, 1e2, 0.001)
soln = odeint(system, initial_state, ktau_vals)

# Obtain solutions and plot
delta = soln[:, 1]
phi = soln[:, 0]

plt.plot(ktau_vals, delta, 'k-')
plt.title(r"$\delta_c$ Behaviour", fontsize=16)
plt.ylabel(r'$\frac{\delta_c}{\mathcal{R}_k}$', fontsize=16, rotation=0, labelpad=10)
plt.xlabel(r'$k\tau$', fontsize=13)
plt.xscale('log')
plt.legend(loc='best')
plt.savefig('delta_behaviour', dpi=300)
plt.show()

plt.plot(ktau_vals, phi, 'k-')
plt.title(r"$\Phi$ Oscillations", fontsize=16)
plt.ylabel(r'$\frac{\Phi}{\mathcal{R}_k}$', rotation=0, fontsize=16)
plt.xlabel(r'$k\tau$', fontsize=13)
plt.xscale('log')
plt.savefig('phi_oscillations', dpi=300)
plt.show()

ref = lambda ktau: 1 / (ktau - 1) ** (1.5)
ref_minus = lambda ktau: -1 / (ktau - 1) ** (1.5)

plt.plot(ktau_vals, [ref(ktau) for ktau in ktau_vals], 'k--')
plt.plot(ktau_vals, [ref_minus(ktau) for ktau in ktau_vals], 'k--')
plt.plot(ktau_vals, phi, 'k-')
plt.ylim(-0.15, 0.15)
plt.xlim(1, 100)
plt.title(r"$\Phi$ Oscillations", fontsize=16)
plt.ylabel(r'$\frac{\Phi}{\mathcal{R}_k}$', rotation=0, fontsize=16)
plt.xlabel(r'$k\tau$', fontsize=13)
plt.xscale('log')
plt.legend(loc='best')
plt.savefig('phi_oscillations_zoomed', dpi=300)
plt.show()
