from q2_b import g_star
import numpy as np
import matplotlib.pyplot as plt

T_rh = 1e14  # GeV
M_pl = 2.4e18

T_a = lambda a: g_star(a) / a

a_t = lambda t: np.sqrt(((np.pi / (3 * M_pl)) * np.sqrt(g_star(t) / 10) * g_star(t) ** (-2 / 3))) * np.sqrt(t)

T_t = lambda t: g_star(a_t(t)) / a_t(t)

a_vals = np.arange(1, 1e6)
t_vals = np.arange(1, 1e6)  # 1/GeV

T_threshold = [100 * 1000, 30 * 1000, 15 * 1000, 1 * 1000, 200, 50, 25]  # MeV
for T_val in T_threshold:
    plt.vlines(T_val, 0, 200, linestyles='dashed', colors='k')

T_a_result = []
for i in range(len(a_vals)):
    T_a_result.append(T_a(a_vals[i]))

a_t_result = []
for i in range(len(t_vals)):
    a_t_result.append(a_t(t_vals[i]))

T_t_result = []
for i in range(len(t_vals)):
    T_t_result.append(T_t(t_vals[i]))

plt.loglog(a_vals, T_a_result)
plt.show()

plt.loglog(a_vals, a_t_result)
plt.show()

plt.loglog(a_vals, T_t_result)
plt.show()
