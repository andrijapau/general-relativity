import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv

ktau_vals = np.arange(1e-2, 1e2, 0.001)

phi = lambda ktau: (3 / 2) * jv(1, ktau) / ktau
ref = lambda ktau: 1 / ktau ** 3

plt.plot(ktau_vals, [phi(ktau) for ktau in ktau_vals], 'b-', label=r'Analytic solution')
plt.plot(ktau_vals, [ref(ktau) for ktau in ktau_vals], 'k--', label=r'Oscillation Envelope')
plt.ylim(-0.5, 2)

plt.ylabel(r'$\frac{\phi}{\mathcal{R}_k}$')
plt.xlabel(r'$k\tau$')
plt.xscale('log')
plt.legend(loc='best')

plt.show()

plt.plot(ktau_vals, [phi(ktau) for ktau in ktau_vals], 'b-', label=r'Analytic solution')
plt.plot(ktau_vals, [ref(ktau) for ktau in ktau_vals], 'k--', label=r'Oscillation Envelope')
plt.ylim(-0.15, 0.15)
plt.xlim(1, 100)
plt.ylabel(r'$\frac{\phi}{\mathcal{R}_k}$')
plt.xlabel(r'$k\tau$')
plt.xscale('log')
plt.legend(loc='best')
plt.show()
