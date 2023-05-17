import time
import numpy as np
import matplotlib.pyplot as plt
from cps import Cps, T_range

def f_cp(cp, T):

    for i, t in enumerate(T_range):
        if i == 0:
            continue
        if T >= T_range[i - 1] and T < T_range[i]:
            i1 = i - 1
            i2 = i
            break

    T1 = T_range[i1]
    T2 = T_range[i2]
    cp1 = cp[i1]
    cp2 = cp[i2]

    return ((T - T1)*cp2 + (T2 - T)*cp1)/(T2 - T1)

def int_cp(cp, Tl, Th, dT=50):
    '''Integrate cp from Tl to Th'''
    dh = 0

    while Tl < Th:
        dh += dT*f_cp(cp, Tl+dT/2)
        Tl += dT
    
    return dh

def find_Tad(cp_p, rhs, tol=100, max_it=50):
    Tmax = 3500
    Tmin = 1000

    T_ad = 2000
    it = 0

    res = rhs - int_cp(cp_p, 298.15, T_ad)

    while abs(res) > tol and it < max_it:
        if res > 0:
            Tmin = T_ad
            T_ad = (Tmin + Tmax)/2
        else:
            Tmax = T_ad
            T_ad = (Tmax + Tmin)/2

        res = rhs - int_cp(cp_p, 298.15, T_ad)
        it += 1

    return T_ad


#       02, H2, CO2, H2O, C2H4, N2
# W =   32, 2, 44, 18, 28, 28 # g/mol

phis = np.arange(0.4, 1.8, 0.2)

W_reac = np.array([2, 32, 28])/1000
W_prod = np.array([18, 2, 32, 28])/1000

Cp_reac = np.array([Cps[1]/W_reac[0], Cps[0]/W_reac[1], Cps[5]/W_reac[2]])
Cp_prod = np.array([Cps[3]/W_prod[0], Cps[1]/W_prod[1], Cps[0]/W_prod[2], Cps[5]/W_prod[3]])

hf_reac = np.array([0, 0, 0]) #J/kg
hf_prod = np.array([-285.83*1000/W_prod[0], 0, 0, 0]) #J/kg

T_ad_hydrogen = np.zeros_like(phis)

for j, phi in enumerate(phis):
    # reactants:
    n_reac = np.array([1, 0.5/phi, 1.88/phi]) # h2, o2, n2
    m_reac = W_reac*n_reac
    m_tot = sum(m_reac)

    #products
    n_prod = np.array([(1 if phi<1 else 1/phi), (0 if phi<1 else 1 - 1/phi), (0.5/phi - 0.5 if phi<1 else 0), 1.88/phi]) # h2o, h2, o2, n2
    m_prod = W_prod*n_prod

    m_frac_reac = m_reac/m_tot
    m_frac_prod = m_prod/m_tot

    Cp_mix_r = np.zeros_like(Cp_reac[0])
    Cp_mix_p = np.zeros_like(Cp_reac[0])

    for i in range(3):
        Cp_mix_r += m_frac_reac[i]*Cp_reac[i]

    for i in range(4):
        Cp_mix_p += m_frac_prod[i]*Cp_prod[i]

    Hr = sum(m_frac_reac*hf_reac) - sum(m_frac_prod*hf_prod)
    Q_r = int_cp(Cp_mix_r, 298.15, 1100)

    T_ad_hydrogen[j] = find_Tad(Cp_mix_p, Q_r+Hr)

# plt.plot(phis, T_ad_hydrogen, marker='o')
# plt.grid()
# plt.xlabel('Equivalence ratio [-]')
# plt.ylabel('T_ad [K]')
# plt.show()
