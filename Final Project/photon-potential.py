import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 16})
V = lambda r, r_s: (1 / r ** 2) * (1 - r_s / r)

r_vals = np.arange(0.1, 5, 0.001)

r_s = 0.5
plt.plot(r_vals, V(r_vals, r_s=r_s), 'k-', linewidth=2, label=r'$V(r)$')
plt.vlines(x=1.5 * r_s, ymin=-5, ymax=2, colors='black', linestyles={"dashed"}, label=r'$r_{\star} = 3M$')
plt.gca().set_aspect('equal')
plt.locator_params(axis='y', nbins=5)
plt.locator_params(axis='x', nbins=5)
plt.ylim(-2, 2)
plt.xlim(0, 4)
plt.legend(loc='best')
plt.title(r"Effective Potential $V(r)$", fontsize=16)
plt.ylabel(r'$V(r)$', rotation=0, fontsize=16, labelpad=10)
plt.xlabel(r'$r$', fontsize=16)
plt.savefig('potential.png', bbox_inches='tight', dpi=400)
plt.show()
