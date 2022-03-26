import matplotlib.ticker
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

# T_evals = arange(0, T_rh)
# plt.loglog(T_evals, g_star(T_evals), 'r-', label=r'$g_{\star}$')
# plt.loglog(T_evals, g_star_s(T_evals), 'b-', label=r'$g_{\star S}$')
# plt.hlines(g_low, min(T_evals), max(T_evals), linestyles='dashed', colors='k', linewidth=0.5)
# plt.hlines(g_high, min(T_evals), max(T_evals), linestyles='dashed', colors='k', linewidth=0.5)
# plt.hlines(g_low_s, min(T_evals), max(T_evals), linestyles='dashed', colors='k', linewidth=0.5)
# plt.title(r'$g(T)$ Transition')
# plt.ylabel(r'$g$')
# plt.xlabel(r'$T$ (GeV)')
# plt.legend(loc='best')
# plt.savefig('g_transition_analytical', dpi=300)
# plt.show()

# T_evals_zoomed = arange(T_0 - 10, T_0 + 100)
# plt.loglog(T_evals_zoomed, g_star(T_evals_zoomed), 'r-')
# plt.loglog(T_evals_zoomed, g_star_s(T_evals_zoomed), 'b-')
# plt.hlines(g_low, min(T_evals_zoomed), max(T_evals_zoomed), linestyles='dashed', colors='k')
# plt.hlines(g_high, min(T_evals_zoomed), max(T_evals_zoomed), linestyles='dashed', colors='k')
# plt.hlines(g_low_s, min(T_evals_zoomed), max(T_evals_zoomed), linestyles='dashed', colors='k')
# plt.show()

# T(a)
a = arange(1, 1e7)


# plt.loglog(a_zoomed, T(a_zoomed), 'k-')
# plt.show()

def T(a):
    return T_rh * (g_high ** (1 / 3)) * (g_star_s(T_rh / a) ** (-1 / 3)) / a


a_crit = T_rh / T_0


def T_approx(a):
    if a < T_rh / T_0:
        return T_rh * (g_high ** (1 / 3)) * (g_high ** (-1 / 3)) / a
    if a > T_rh / T_0:
        return T_rh * (g_high ** (1 / 3)) * (g_low_s ** (-1 / 3)) / a


# T_approx_array = [T_approx(a_val) for a_val in a]
# T_analytic = [T(a_val) for a_val in a]
# plt.loglog(a, T_analytic, 'k-', linewidth=0.5, label="Analytic")
# plt.loglog(a, T_approx_array, 'b-', linewidth=1, label='Piecewise Approx.')
# plt.title(r'$T(a)$')
# plt.ylabel(r'$T$ GeV')
# plt.xlabel(r'$a$')
# plt.legend(loc='best')
# plt.savefig('temperature', dpi=300)
# plt.show()
# #
# plt.loglog(a, T_analytic, 'k-', linewidth=0.5, label="Analytic")
# plt.loglog(a, T_approx_array, 'b-', linewidth=1, label='Piecewise Approx.')
# plt.title(r'$T(a)$')
# plt.ylabel(r'$T$ GeV')
# plt.xlabel(r'$a$')
# plt.legend(loc='best')
# plt.xlim(a_crit - 900, a_crit + 9000)
# plt.ylim(T_0 - 90, T_0 + 900)
# plt.savefig('temperature_zoomed', dpi=300)
# plt.show()

print("Ratio at a=1e7: ", T_approx(1e7) / T(1e7))

# # Solve FRW equation
# t_rh = ((3 * M_pl * sqrt(10)) / (2 * pi)) * 1 / (T_rh ** 2 * sqrt(g_high))
# t_rh = (1 / T_rh) / (1.5e24)
#
# t_evals = arange(t_rh, 1e12 * t_rh, step=1e6 * t_rh)
# a_dot = lambda t, a: (((T_rh ** 2 * pi) / (3 * M_pl * sqrt(10))) * (
#         g_high ** (2 / 3) * g_star((1 / (t * 1.5e24))) ** (1 / 2)) / (
#                               g_star_s((1 / (t * 1.5e24))) ** (2 / 3))) * (1 / a)
# soln = odeint(a_dot, y0=1, t=t_evals, tfirst=True)
# print(soln)
#
# plt.loglog(t_evals, soln)
# plt.show()


# t_rh = ((3 * M_pl * sqrt(10)) / (2 * pi)) * 1 / (T_rh ** 2 * sqrt(g_high))
t_rh = (1 / T_rh) / (1.5e24)
print(t_rh)
t_evals = arange(t_rh, 1e12 * t_rh, step=1e5 * t_rh)

t0 = (1 / T_0) / (1.5e24)
print(t0)


def a_piecewise(t):
    C_high = ((T_rh ** 2 * pi) / (3 * M_pl * sqrt(10))) * sqrt(g_high)
    C_low = ((T_rh ** 2 * pi) / (3 * M_pl * sqrt(10))) * (g_high ** (2 / 3) * sqrt(g_low)) / (g_low_s ** (2 / 3))
    if t < t0:
        # return sqrt(2 * C_high) * sqrt(t)
        return sqrt(t / t_rh)
    if t > t0:
        # return sqrt(2 * C_low) * sqrt(t)
        return ((g_high ** (1 / 6) * g_low ** (1 / 2)) / (g_low_s ** (2 / 3))) * sqrt(t / t_rh)


a_soln = [a_piecewise(t_val) for t_val in t_evals]
# plt.loglog(t_evals, a_soln, 'k-', label=r'Approximate', linewidth=1.5)
# plt.vlines(t0, min(a_soln), max(a_soln), linestyles='dashed', colors='k', linewidth=0.5)
# plt.title(r'$a(t)$ Transition')
# plt.ylabel(r'$a$')
# plt.xlabel(r'$t$ (s)')
# plt.legend(loc='best')
# plt.savefig('a_transition_piecewise', dpi=300)
# plt.show()
# T_t_soln = [T_approx(a_val) for a_val in a_soln]
# plt.loglog(t_evals, T_t_soln, 'k-', label=r'Approximate', linewidth=1.5)
# plt.title(r'$T(t)$ Transition')
# plt.ylabel(r'$T$')
# plt.xlabel(r'$t$ (s)')
# plt.legend(loc='best')
# plt.savefig('T_t_function', dpi=300)
# plt.show()
