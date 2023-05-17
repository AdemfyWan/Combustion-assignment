# combustion assignment 1, part 3
# Minrui Wan, Georgi


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline, make_interp_spline
from scipy import optimize
from scipy.optimize import fsolve


def part3():


    df_C2H4_1 = pd.read_csv('C2H4_part3.csv',index_col=0)
    df_C2H4_mat1 = df_C2H4_1.to_numpy()

    df_C2H4_10 = pd.read_csv('C2H4_part3_10bar.csv', index_col=0)
    df_C2H4_mat10 = df_C2H4_10.to_numpy()
    # Numerator
    NO_cea_1  = []
    temp_1= []
    phi_1 = []
    NO_cea_10 = []
    temp_10 = []
    phi_10 = []
    for i in range(4):
        NO_cea_1.append(df_C2H4_mat1[i][13])
        temp_1.append(df_C2H4_mat1[i][0])
        phi_1.append(df_C2H4_mat1[i][1])
        NO_cea_10.append(df_C2H4_mat10[i][13])
        temp_10.append(df_C2H4_mat10[i][0])
        phi_10.append(df_C2H4_mat10[i][1])

    def func1_1(z):
        phi = 0.4
        Kp = 3.0578446650349116e-05
        return 4 * z ** 2 - Kp * (11.28 / phi - z) * (3 / phi - 3 - z)
    def func1_10(z):
        phi = 0.4
        Kp = 3.066419019332548e-05
        return 4 * z ** 2 - Kp * (11.28 / phi - z) * (3 / phi - 3 - z)

    def func2_1(z):
        phi = 0.6
        Kp = 0.0004406411875594603
        return 4 * z ** 2 - Kp * (11.28 / phi - z) * (3 / phi - 3 - z)
    def func2_10(z):
        phi = 0.6
        Kp = 0.0004565957562772559
        return 4 * z ** 2 - Kp * (11.28 / phi - z) * (3 / phi - 3 - z)

    def func3_1(z):
        phi = 0.8
        Kp = 0.0017742244560127788
        return 4 * z ** 2 - Kp * (11.28 / phi - z) * (3 / phi - 3 - z)
    def func3_10(z):
        phi = 0.8
        Kp = 0.0021215611103635947
        return 4 * z ** 2 - Kp * (11.28 / phi - z) * (3 / phi - 3 - z)

    def func4_1(z):
        phi = 0.99
        Kp = 0.0031294461457621776
        return 4 * z ** 2 - Kp * (11.28 / phi - z) * (3 / phi - 3 - z)
    def func4_10(z):
        phi = 0.99
        Kp = 0.004341702369050261
        return 4 * z ** 2 - Kp * (11.28 / phi - z) * (3 / phi - 3 - z)
    # Combustion equation
    #[C2H4] + 3/phi*([O2] + 3.76*[N2]) ==== 2[CO2] +  2[H2O] + 3/phi*3.76[N2] + (3/phi-3)[O2]
    # [N2]  + [O2]  ==== 2[NO]
    # gibbs energy
    Hf_NO   =  90.29  #kJ/mol
    Sf_NO   =  210.66 #J/mol/K
    Sf_O2   = 205.04  #J/mol/K
    Sf_N2   = 191.50  #J/mol/K
    R0 = 8.3145*10**(-3)  # [J/mol/K]
    dG_R0_1 = []
    dG_R0_10= []
    dH = []
    dS = []
    NO_eq_1 = []
    NO_eq_10= []
    Kp_1   = []
    Kp_10  = []
    n_tot= []
    for i in range(4):
        dH.append(2 * Hf_NO)
        dS.append((2 * Sf_NO - Sf_N2 - Sf_O2) / 1000)
        dG_R0 = ((2 * Hf_NO) - temp_1[i] * (2 * Sf_NO - Sf_N2 - Sf_O2)/1000)
        dG_R0_1.append(dG_R0)
        Kp_1.append(np.exp(-dG_R0_1[i] / R0 / temp_1[i]))
        n_tot.append(4 + 3 * 3.76/phi_1[i]+3/phi_1[i]-3)
        #z = bisect(f,0.00000001,1)
        if i == 0:
            z = fsolve(func1_1,0)
        elif i == 1:
            z = fsolve(func2_1, 0)
        elif i == 2:
            z = fsolve(func3_1, 0)
        elif i == 3:
            z = fsolve(func4_1, 0)
        NO_eq_1.append(2*z[0]/n_tot[i])
    for i in range(4):
        dG_R0 = ((2 * Hf_NO) - temp_10[i] * (2 * Sf_NO - Sf_N2 - Sf_O2)/1000)
        dG_R0_10.append(dG_R0)
        Kp_10.append(np.exp(-dG_R0_10[i] / R0 / temp_10[i]))
        n_tot.append(4 + 3 * 3.76/phi_1[i]+3/phi_1[i]-3)
        #z = bisect(f,0.00000001,1)
        if i == 0:
            z = fsolve(func1_10,0)
        elif i == 1:
            z = fsolve(func2_10, 0)
        elif i == 2:
            z = fsolve(func3_10, 0)
        elif i == 3:
            z = fsolve(func4_10, 0)
        NO_eq_10.append(2*z[0]/n_tot[i])

    print('dH:\n',dH)
    print('dS:\n',dS)
    print('ntot:\n',n_tot)
    print('dG_R0_1bar:\n',dG_R0_1)
    print('dG_R0_10bar:\n',dG_R0_10)
    print('temp_1bar:\n',temp_1)
    print('temp_10bar:\n',temp_10)
    print('Kp_1bar:\n',Kp_1)
    print('Kp_10bar:\n',Kp_10)
    print()
    print('equilibrium result @ 1bar:\n',NO_eq_1)
    print('equilibrium result @ 10bar:\n',NO_eq_10)
    print()
    print('CEA result @ 1bar:\n',NO_cea_1)
    print('CEA result @ 10bar:\n',NO_cea_10)
#----------------------------------------------------------------
    plotboth = 1

    AxisX = [0.4, 0.6, 0.8, 0.99]
    AxisXnew = np.linspace(0.4, 0.99, 150)
    spline_eq1 = make_interp_spline(AxisX, NO_eq_1, k=3)
    spline_CEA1 = make_interp_spline(AxisX, NO_cea_1, k=3)
    NO_eq_smth1 = spline_eq1(AxisXnew)
    NO_cea_smth1 = spline_CEA1(AxisXnew)
    spline_eq10 = make_interp_spline(AxisX, NO_eq_10, k=3)
    spline_CEA10 = make_interp_spline(AxisX, NO_cea_10, k=3)
    NO_eq_smth10 = spline_eq10(AxisXnew)
    NO_cea_smth10 = spline_CEA10(AxisXnew)
    plt.figure('NO mass fraction (1bar, 600K)')
    plt.plot(AxisX, NO_eq_1, 'o',color = 'k')
    plt.plot(AxisX, NO_cea_1, 'o',color = 'g')
    plt.plot(AxisXnew, NO_eq_smth1, label='eq. @ 1bar', color='k')
    plt.plot(AxisXnew, NO_cea_smth1, label='cea @ 1bar', color='g')
    if plotboth == 1:
        plt.plot(AxisXnew, NO_eq_smth10, label='eq. @ 10bar', color='k', ls= '--')
        plt.plot(AxisXnew, NO_cea_smth10, label='cea @ 10bar', color='g', ls= '--')
    # plt.plot(AxisX, NO_H2, label='NO_H2')
    plt.legend()
    plt.xlabel('Equivalence Ratio')
    plt.ylabel('NO mass fraction')
    plt.show()


part3()