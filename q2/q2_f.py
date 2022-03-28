# Import modules
from q2_c import *


def s(T):
    # Calculate entropy at a given temperature
    return ((2 * np.pi ** 2) / 45) * np.sqrt(g_star_s(T)) * (T) ** 3


def calculate_gw(omega, rho, m, s):
    # Calculates gw given various parameters
    return ((6e-17) / ((omega * rho) / (m * s))) ** (1 / 3.8)


# Find stuff at ep annihilation temperature
T_ep = 2e-7  # GeV
H0_GeV = (np.pi / (3 * np.sqrt(10))) * np.sqrt(g_star(T_ep)) * (T_ep ** 2) / M_pl  # GeV
inverse_gev_to_seconds = (1.52e24)
H0_seconds = (1 / H0_GeV) / inverse_gev_to_seconds
G = 6.67e-8  # grams/ cm3 s2
M_pl_sq = 1 / (8 * np.pi * G)

# Parameters needed to calculate gw
gcm3_to_gev4 = 2.32e17
critical_rho_gev4 = (1 / ((3 * H0_seconds ** 2) / M_pl_sq)) / gcm3_to_gev4
omega_x = 0.01
s_today = s(T_ep)
m_x = 500

# Print result
print("gw: {}".format(round(calculate_gw(omega_x, critical_rho_gev4, m_x, s_today), 3)))
