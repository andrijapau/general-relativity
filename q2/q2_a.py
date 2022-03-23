from numpy import arctan, pi, arange, sqrt
import matplotlib.pyplot as plt
from scipy.integrate import odeint

T_rh = 1e5
g_high = 500
g_low = 50
g_low_s = 200
T_0 = 100
T_width = 3
M_pl = 2.4e18

# Transition in g(T)

g_star = lambda T: g_low + ((g_high - g_low) / pi) * (arctan((T - T_0) / T_width) + pi / 2)
g_star_s = lambda T: g_low_s + ((g_high - g_low_s) / pi) * (arctan((T - T_0) / T_width) + pi / 2)
#
# T_evals = arange(0, T_rh)
# plt.loglog(T_evals, g_star(T_evals), 'r-')
# plt.loglog(T_evals, g_star_s(T_evals), 'b-')
# plt.hlines(g_low, min(T_evals), max(T_evals), linestyles='dashed', colors='k')
# plt.hlines(g_high, min(T_evals), max(T_evals), linestyles='dashed', colors='k')
# plt.hlines(g_low_s, min(T_evals), max(T_evals), linestyles='dashed', colors='k')
# plt.show()
#
# T_evals_zoomed = arange(T_0 - 10, T_0 + 100)
# plt.loglog(T_evals_zoomed, g_star(T_evals_zoomed), 'r-')
# plt.loglog(T_evals_zoomed, g_star_s(T_evals_zoomed), 'b-')
# plt.hlines(g_low, min(T_evals_zoomed), max(T_evals_zoomed), linestyles='dashed', colors='k')
# plt.hlines(g_high, min(T_evals_zoomed), max(T_evals_zoomed), linestyles='dashed', colors='k')
# plt.hlines(g_low_s, min(T_evals_zoomed), max(T_evals_zoomed), linestyles='dashed', colors='k')
# plt.show()

# T(a)
a = arange(1, 1e7)
a_zoomed = arange(T_0 - 50, T_0 + 50)

T = lambda a: (g_star_s(a) ** (-1 / 3)) / a

# plt.loglog(a_zoomed, T(a_zoomed), 'k-')
# plt.show()
plt.loglog(a, T(a), 'k-')


def T_approx(a):
    if a > T_0:
        return (g_high ** (-1 / 3)) / a
    if a < T_0:
        return (g_low_s ** (-1 / 3)) / a


T_approx_array = [T_approx(a_val) for a_val in a]
T_approx_array_zoomed = [T_approx(a_val) for a_val in a_zoomed]

plt.loglog(a, T_approx_array, 'b-')
plt.show()

plt.loglog(a_zoomed, T_approx_array_zoomed, 'b-')
plt.loglog(a_zoomed, T(a_zoomed), 'k-')
plt.show()

print("Ratio at a=1e7: ", T_approx(1e7) / T(1e7))

# Solve FRW equation
t_rh = 1

t_evals = arange(t_rh, 1e5 * t_rh)
a_dot = lambda t, a: (pi / 3) * sqrt(g_star(a) / 10) * (1 / M_pl) * T(a) ** 2
soln = odeint(a_dot, y0=1, t=t_evals, tfirst=True)

plt.loglog(t_evals, soln)
plt.show()
