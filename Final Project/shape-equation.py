import numpy as np
from sympy import atan2, acos
import matplotlib.pyplot as plt
import matplotlib

plt.rcParams.update({'font.size': 16})

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
from scipy.integrate import odeint, RK45


def dU_dx(U, x):
    if (U[0] < 0.0001):
        return [0, 0]
    return [U[1], 1.5 * U[0] ** 2 - U[0]]


def F(y, u, x):
    return 1.5 * y ** 2 - y


def asSpherical(x, y, z):
    # takes list xyz (single coord)
    r = np.sqrt(x * x + y * y + z * z)
    theta = acos(z / r)
    phi = atan2(y, x)
    return r, theta, phi


M = 1

b = 1
a = 0

r_s = 2 * M
delta0 = (0) * np.pi / 180
r0, phi0, theta0 = 10 * r_s, 0, np.pi / 2 - delta0

R = 5 * M
U0 = [1 / r0, 0]

a = 0
b = np.pi - delta0
N = 1000
xpoints = np.linspace(a, b, N)

h = (b - a) / N

ypoints = []
upoints = []

y = 1 / r0
u = 0

for x in xpoints:
    if y < 100:
        y = 0.001
        u = 0.001

    ypoints.append(y)
    upoints.append(u)

    m1 = h * u
    k1 = h * F(y, u, x)  # (x, v, t)

    m2 = h * (u + 0.5 * k1)
    k2 = h * F(y + 0.5 * m1, u + 0.5 * k1, x + 0.5 * h)

    m3 = h * (u + 0.5 * k2)
    k3 = h * F(y + 0.5 * m2, u + 0.5 * k2, x + 0.5 * h)

    m4 = h * (u + k3)
    k4 = h * F(y + m3, u + k3, x + h)

    y += (m1 + 2 * m2 + 2 * m3 + m4) / 6
    u += (k1 + 2 * k2 + 2 * k3 + k4) / 6

r = [1 / y for y in ypoints]
plt.plot(xpoints, ypoints)
plt.show()

# Us = odeint(dU_dx, U0, xs)
# ys = Us[:, 0]
plt.plot(xpoints, r, 'k-', linewidth=2)
plt.title(r"Shape Equation", fontsize=16)
plt.ylabel(r'$u(\varphi)$', rotation=0, fontsize=16, labelpad=10)
plt.xlabel(r'$\varphi$', fontsize=16)
plt.savefig('shape-equation.png', bbox_inches='tight', dpi=400)
plt.ylim(0, r0)
plt.show()

fig = plt.figure()
ax = plt.axes(projection="3d")

x = r * np.sin(xpoints)
y = [0] * len(r)
z = r * np.cos(xpoints)
ax.plot3D(x, y, z, 'gray')
ax.scatter3D(r0 * np.sin(theta0) * np.cos(phi0), r0 * np.sin(theta0) * np.sin(phi0), r0 * np.cos(theta0))
ax.scatter3D(0, 0, 0)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
