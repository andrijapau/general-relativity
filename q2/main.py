import numpy as np
import matplotlib.pyplot as plt

# (iv)

T_rh = 10e5  # GeV
g_high = 500
g_low = 50
g_low_s = 200
T_0 = 100  # GeV
T_width = 3  # GeV
M_pl = 2.4e18  # GeV


def g(T, T_0, T_width, g_low, g_high):
    return g_low + (g_high - g_low / np.pi) * (np.arctan((T - T_0) / T_width) + np.pi / 2)
