import numpy as np
import matplotlib.pyplot as plt

R = 8.314

def cp(ah, al, T, T_mid):
    if T > T_mid:
        a = ah
    else:
        a = al
    return (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4)*R

get_cp = np.vectorize(cp, excluded=['ah', 'al'])

with open('species.dat', 'r') as file:
    data = file.readlines()

names = []
Ts = []
a_high = []
a_low = []
cps = []


for i, row in enumerate(data):
    els = row.split()
    if i%3 == 0:
        names.append(els[0])
        Ts.append([float(num) for num in els[1:]])
    elif i%3 == 1:
        a_high.append([float(num) for num in els])
    else:
        a_low.append([float(num) for num in els])

T_range = np.linspace(290, 3600)

for i, name in enumerate(names):
    cps.append(get_cp(ah=a_high[i], al=a_low[i], T=T_range, T_mid=Ts[i][2]))

Cps = np.array(cps) # J/mol K
W = np.array([32, 2, 44, 18, 28, 28])/1000



if __name__ == '__main__':

    for i, name in enumerate(names):
        Cp = Cps[i]/W[i]
        plt.plot(T_range, Cp, label=name)


    plt.legend()
    plt.ylabel('Cp[J/kg K]')
    plt.xlabel('T[K]')
    plt.grid()
    plt.show()
