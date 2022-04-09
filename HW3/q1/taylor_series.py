from sympy import symbols, exp, diff, simplify, Function, integrate, sqrt, pi

# A = 45 / (4 * pi ** 4)
# g_x = 4
# g_s = 106.75

X, N, l, N_eq, dN, dX, a, A = symbols('X N \u03BB N_eq dN dX a A')
# N_eq = Function('N_eq')(X)
# f_p = a ** 2 / (exp(sqrt(a ** 2 + exp(X) ** 2)) - 1)
# N_eq = A * integrate(f_p, a)
# dN_eq_dX = diff(N_eq, X)
# print(dN_eq_dX)
F = -(l / exp(X)) * (exp(N) - exp(2 * N_eq - N))
dF_dN = diff(F, N)
dF_dX = diff(F, X)
print(dF_dN)
print(dF_dX)
F_approx = simplify(F + dF_dN * dN + dF_dX * dX)
print(F_approx)
