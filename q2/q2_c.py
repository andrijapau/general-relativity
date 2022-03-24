from q2_b import *

T_rh = 1e14 # GeV

T_a = lambda a: g_star(a) / a 

a_t = lambda t: sqrt(((pi/(3*M_pl)) * sqrt(g_star(t)/10) * g_star(t)**(-2/3))) * sqrt(t) 

T_t = lambda t: g_star(a(t)) / a(t)

a_vals = np.arange(1e1, 1e9, 10) 
#t_vals = np.arange(1e12, 1e26) # 1/GeV


plt.loglog(a_vals, T_a(a_vals))
plt.show()
