from q2_c import *

# ============================================
# Find T(t) from t_rh to today!
# ============================================

logt_today = np.log(6.62e41)
T_today = T_approx_exp(a_piecewise(logt_today))

# ============================================
# Find H0
# ============================================

H0_GeV = (np.pi / (3 * np.sqrt(10))) * np.sqrt(g_star(T_today)) * (T_today ** 2) / M_pl  # GeV

print("H0 (GeV): {}".format(H0_GeV))

# ============================================
# Find critical density for universe
# ============================================

# Conversion Factor
inverse_gev_to_seconds = (1.52e24)

# Convert H0 to seconds
H0_seconds = (1 / H0_GeV) / inverse_gev_to_seconds
print("H0 (Seconds): {}".format(H0_seconds))

# Convert M_pl to si units
G = 6.67e-8  # grams/ cm3 s2
M_pl_sq = 1 / (8 * np.pi * G)

# Calculate critical density
critical_rho = 1 / ((3 * H0_seconds ** 2) / M_pl_sq)

print("Critical Density (g / cm3): {}".format(critical_rho))

# ============================================
# Calculate Vacuum Density
# ============================================

omega_vac = 0.68
rho_DM_gcm3 = omega_vac * critical_rho

print("Vacuum Density (g / cm3): {}".format(rho_DM_gcm3))

gcm3_to_ev4 = 2.32e17 * (1 / (1e9) ** 4)
rho_DM_ev4 = rho_DM_gcm3 / gcm3_to_ev4
print("Vacuum Density (eV^4): {}".format(rho_DM_ev4))

# ============================================
# Calculate Cosmological Constant Today
# ============================================

c = 1e8
CC = rho_DM_gcm3 * (8 * np.pi * G / c ** 2)
print("Cosmological Constant Today (1/m^2): {}".format(CC))

# ============================================
# Calculate Cosmological Constant before EW
# ============================================

# Same process as before but with different T0!!
T0_EW = 100  # GeV
H0_GeV_EW = (np.pi / (3 * np.sqrt(10))) * np.sqrt(g_star(T0_EW)) * (T0_EW ** 2) / M_pl  # GeV
inverse_gev_to_seconds = (1.52e24)
H0_seconds = (1 / H0_GeV_EW) / inverse_gev_to_seconds
G = 6.67e-8  # grams/ cm3 s2
M_pl_sq = 1 / (8 * np.pi * G)
critical_rho = 1 / ((3 * H0_seconds ** 2) / M_pl_sq)
omega_vac = 0.68
rho_DM_gcm3 = omega_vac * critical_rho
c = 1e8
CC_EW = rho_DM_gcm3 * (8 * np.pi * G / c ** 2)

print("Ratio b/w CC's: {}".format(CC_EW / CC))
