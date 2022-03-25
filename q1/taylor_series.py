from sympy import symbols, exp, diff, simplify

X, N, l, N_eq, dN, dX = symbols('X N \u03BB N_eq dN dX')

F = -(l / exp(X)) * (exp(N) - exp(2 * N_eq - N))
dF_dN = diff(F, N)
dF_dX = diff(F, X)
print(dF_dN)
print(dF_dX)
F_approx = simplify(F + dF_dN * dN + dF_dX * dX)
print(F_approx)
