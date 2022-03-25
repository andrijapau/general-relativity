from q2_b import g_star, g_star_s
import numpy as np
import matplotlib.pyplot as plt

T_rh = 1e14  # GeV
M_pl = 2.4e18
one_second = 1.52e24  # 1/Gev

T_a = lambda a: g_star_s(a) ** (-1 / 3) / a

a_t = lambda t: np.sqrt(((np.pi / (3 * M_pl)) * np.sqrt(g_star(t) / 10) * g_star_s(t) ** (-2 / 3))) * np.sqrt(t)

T_t = lambda t: g_star(a_t(t)) / a_t(t)

a_vals = np.arange(1, 1e7)
t_vals = np.arange(1e-12, 1e2)  # 1/GeV

T_threshold = [100 * 1000, 30 * 1000, 15 * 1000, 1 * 1000, 200, 50, 25]  # MeV
for T_val in T_threshold:
    plt.vlines(T_val, 0, 200, linestyles='dashed', colors='k', linewidth=0.75)

T_a_result = []
for i in range(len(a_vals)):
    T_a_result.append(T_a(a_vals[i]))

a_t_result = []
for i in range(len(t_vals)):
    a_t_result.append(a_t(t_vals[i]))

T_t_result = []
for i in range(len(t_vals)):
    T_t_result.append(T_t(t_vals[i]))

plt.loglog(a_vals, T_a_result, 'k-', linewidth=1.5)
plt.title(r'$T(a)$')
plt.xlabel(r'$a$')
plt.ylabel(r'$T(a)$')
plt.axis([min(a_vals), max(a_vals), min(T_a_result), max(T_a_result)])
plt.show()

plt.loglog(t_vals, a_t_result)
plt.show()

plt.loglog(t_vals, T_t_result)
plt.show()
