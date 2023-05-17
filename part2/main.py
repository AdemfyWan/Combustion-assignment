# combustion assignment 1, part 2
# Minrui Wan, Georgi


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline, make_interp_spline
from openpyxl import load_workbook

def main():
    #wb = load_workbook(filename='cea-H2.xlsx')
    #ws = wb['H2_600']
    #all_rows = list(ws.rows)
    #print(all_rows[:5])

    df_H2 = pd.read_csv('H2_600.csv', index_col=0)
    df_C2H4=pd.read_csv('C2H4_600.csv',index_col=0)
    # for p = 1 bar
    for j in range(3):
        df_H2_t = df_H2.iloc[j::3,:]
        df_H2_mat = df_H2_t.to_numpy()
        temp = []
        for i in range(8):
            temp.append(df_H2_mat[i][0])
        if j == 0:
            temp_p1 = temp
        elif j == 1:
            temp_p10 = temp
        else:
            temp_p40 = temp


    AxisX = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8]
    AxisXnew = np.linspace(0.4, 1.8, 300)
    spline_p1 = make_interp_spline(AxisX, temp_p1,k=3)
    spline_p10= make_interp_spline(AxisX, temp_p10,k=3)
    spline_p40= make_interp_spline(AxisX, temp_p40,k=3)
    T_p1_smth = spline_p1(AxisXnew)
    T_p10_smth= spline_p10(AxisXnew)
    T_p40_smth= spline_p40(AxisXnew)
    plt.figure('Variation of Adiabatic flame temperature at different pressure and ER')
    plt.plot(AxisX, temp_p1,'o')
    plt.plot(AxisX, temp_p10,'o')
    plt.plot(AxisX, temp_p40, 'o')
    plt.plot(AxisXnew,T_p1_smth,label='p = 1bar')
    plt.plot(AxisXnew,T_p10_smth,label='p = 10bar')
    plt.plot(AxisXnew, T_p40_smth,label='p = 40bar')
    plt.legend()
    plt.xlabel('equivalence ratio')
    plt.ylabel('Temperature')
    plt.show()



    # NOx emission
    df_H2_p1 = df_H2.iloc[::3,:]
    df_H2_p1_mat = df_H2_p1.to_numpy()
    NO_H2 = []
    NO2_H2= []
    NOx_H2= []
    for i in range(8):
        NO_H2.append(df_H2_p1_mat[i][13])
        NO2_H2.append(df_H2_p1_mat[i][14])
        NOx_H2.append(NO_H2[i] + NO2_H2[i])
    df_C2H4_p1 = df_C2H4.iloc[::3,:]
    df_C2H4_p1_mat = df_C2H4_p1.to_numpy()
    NO_C2H4 = []
    NO2_C2H4 = []
    NOx_C2H4 = []
    for i in range(8):
        NO_C2H4.append(df_C2H4_p1_mat[i][13])
        NO2_C2H4.append(df_C2H4_p1_mat[i][14])
        NOx_C2H4.append(NO_C2H4[i] + NO2_C2H4[i])
    spline_Nox_H2   = make_interp_spline(AxisX,NOx_H2, k=3)
    spline_NOx_C2H4 = make_interp_spline(AxisX,NOx_C2H4, k=3)
    NOx_C2H4_smth = spline_NOx_C2H4(AxisXnew)
    NOx_H2_smth   = spline_Nox_H2(AxisXnew)
    plt.figure('NOx emission (1bar, 600K)')
    plt.plot(AxisX, NOx_C2H4,'o')
    plt.plot(AxisX, NOx_H2,'o')
    plt.plot(AxisXnew, NOx_C2H4_smth, label='NOx_C2H4')
    #plt.plot(AxisX, NO_C2H4, label='NO_C2H4')
    plt.plot(AxisXnew, NOx_H2_smth, label='NOx_H2')
    #plt.plot(AxisX, NO_H2, label='NO_H2')
    plt.legend()
    plt.xlabel('equivalence ratio')
    plt.ylabel('N emission')
    plt.show()

    print(NO_C2H4)

main()