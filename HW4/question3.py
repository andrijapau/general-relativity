# Import modules
import numpy as np
from scipy.integrate import odeint


def a_dot(t, a):
    '''Solves a_dot in coordinate time'''
    omega_r0 = 9.4e-5
    omega_m0 = 0.32
    omega_c0 = 1 - omega_m0 - omega_r0
    omega_k0 = 0
    a0 = 1
    h = 0.67
    numsecondsinGy = 1e9 * 3.154e7
    H0 = (h / 3.086e17) * numsecondsinGy
    return H0 * a * np.sqrt((omega_r0 * (a0 / a) ** 4 + omega_m0 * (a0 / a) ** 3 + omega_k0 * (
            a0 / a) ** 3 + omega_c0))


def H(a):
    '''Solves for H at a given a using FRW equation'''
    omega_r0 = 9.4e-5
    omega_m0 = 0.32
    omega_c0 = 1 - omega_m0 - omega_r0
    omega_k0 = 0
    a0 = 1
    h = 0.67
    c = 3e5  # Mpc / s
    H0 = (100 * h) / c
    return H0 * np.sqrt((omega_r0 * (a0 / a) ** 4 + omega_m0 * (a0 / a) ** 3 + omega_k0 * (
            a0 / a) ** 3 + omega_c0))


# Define t space
t_evals = np.linspace(0.001, -13.81964, 1000000)

# Solve FRW
soln_backward = odeint(a_dot, y0=[1, 1], t=t_evals, tfirst=True)

# Solve for conformal time
dt = t_evals[1] - t_evals[0]
a_inverse_int = np.trapz(1 / soln_backward, dx=-dt, axis=0)
print("Age of Universe in Conformal Time: {:f}e10 Years".format(a_inverse_int[0] / 10))

# Solve for a_eq
matter_radiation_time = 0.00006  # GYr
t_evals = np.linspace(0.001, -13.81965 + matter_radiation_time, 1000000)
dt = t_evals[1] - t_evals[0]
soln_backward = odeint(a_dot, y0=[1, 1], t=t_evals, tfirst=True)
print("a at MR equality: {:2e}".format(soln_backward[-1][-1]))

# Solve for H at various a values
print("H(a_eq) = {}".format(H(2.94e-4)))
print("H(a = 0.5) = {}".format(H(0.5)))
