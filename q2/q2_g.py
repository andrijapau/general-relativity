# Import modules
from q2_c import *


def s(T):
    # Calculate entropy at a given temperature
    return ((2 * np.pi ** 2) / 45) * np.sqrt(g_star_s(T)) * (T) ** 3


def calculate_omega(g_w, rho, m, s):
    # Calculates Omega given various parameters
    return (m * s * (6e-17) * g_w ** (-3.8)) / (rho)


# Calculate critical rho today
today_in_seconds = 6.62e41
logt_today = np.log(6.62e41)
T_today = T_approx_exp(a_piecewise(logt_today))
H0_GeV = (np.pi / (3 * np.sqrt(10))) * np.sqrt(g_star(T_today)) * (T_today ** 2) / M_pl  # GeV
inverse_gev_to_seconds = (1.52e24)
H0_seconds = (1 / H0_GeV) / inverse_gev_to_seconds
G = 6.67e-8  # grams/ cm3 s2
M_pl_sq = 1 / (8 * np.pi * G)

# Parameters needed to calculate gw
gcm3_to_gev4 = 2.32e17
critical_rho_gev4 = (1 / ((3 * H0_seconds ** 2) / M_pl_sq)) / gcm3_to_gev4
g_w = 0.089
s_today = s(T_today)
m_x = 500

# Print result
print("Omega: {}".format(round(calculate_omega(g_w, critical_rho_gev4, m_x, s_today), 3)))
