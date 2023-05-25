# 5.24
# author: Minrui Wan
import numpy as np
import matplotlib.pyplot as plt

#constants
P = 101325 #[Pa]
T = 300 #[K]
R = 8.3145  #[J/K/mol]
W_H2 = 2  #[g/mol]
W_O2 = 32 #[g/mol]
W_N2 = 28 #[g/mol]
lambda_H2 = 0.186 #[mW/m/K]
lambda_O2 = 0.02674 #[]
lambda_N2 = 0.02597 #[]
M_HO = 7.88 * 10**(-5) #[m^2/s]
M_HN = 7.47 * 10**(-5) #[]
M_NO = 2.06 * 10**(-5)
M_OO = 2.07 * 10**(-5)
M_NN = 2.04 * 10**(-5)
bound = 0.4/14
Cp_H2  = 14.30  #[J/g/K]
Cp_O2  = 0.918  #[J/g/K]
Cp_N2  = 1.040  #[]

# mole fraction
def X_O2(x):
    return 0.4 - 4*x

def X_H2(x):
    return 0.4 - 4*x

def X_N2(x):
    return 10*x

# molar mass of mixture
def W_k(x):
    return X_O2(x)* W_H2 + X_O2(x)* W_O2 + X_N2(x)* W_N2

# mass fraction
def Y_k_x1(x):
    return X_O2(x) * W_H2 / W_k(x)

def Y_k_x2(x):
    return X_O2(x) * W_O2 / W_k(x)

def Y_k_x3(x):
    return X_N2(x) * W_N2 / W_k(x)

# density
def rho(x):
    return (P * W_k(x)) / (R * T)

# thermal conductivity of mixture
def lamda(x):
    return 1/2* ((X_O2(x)*0.186 + X_O2(x)*0.02674 + X_N2(x)*0.02597)+
                 (X_O2(x)/0.186 + X_O2(x)/0.02674 + X_N2(x)/0.02597)**(-1))

# Model 1:
def model1_O2(x):
    if x < bound:
        return M_OO
    else:
        return M_NO

def model1_H2(x):
    if x < bound:
        return M_HO
    else:
        return M_HN

def model1_N2(x):
    if x < bound:
        return M_NO
    else:
        return M_NN

# Model 2:
def model2_H2(x):
    return 1/((X_O2(x)/(1-X_H2(x))) / M_HO + (X_N2(x)/(1-X_H2(x))) / M_HN)

def model2_O2(x):
    return 1/((X_H2(x)/(1-X_O2(x))) / M_HO + (X_N2(x)/(1-X_O2(x))) / M_NO)

def model2_N2(x):
    return 1/((X_H2(x)/(1-X_N2(x))) / M_HN + (X_O2(x)/(1-X_N2(x))) / M_NO)

# Model 3:
def model3_H2(x):
    return lamda(x) / rho(x) / Cp_H2* 10**(-6)

def model3_O2(x):
    return lamda(x) / rho(x) / Cp_O2* 10**(-6)

def model3_N2(x):
    return lamda(x) / rho(x) / Cp_N2* 10**(-6)

# Model 4:
def model4_H2(x):
    return lamda(x) / 0.30 / rho(x) / Cp_H2 * 10**(-6)

def model4_O2(x):
    return lamda(x) / 1.11 / rho(x) / Cp_O2 * 10**(-6)

def model4_N2(x):
    return lamda(x) / 1.00/ rho(x) / Cp_N2 * 10**(-6)

def part1():

    x = np.linspace(0, 0.1, 150)
    #X_O2 = 0.4 - 4*x
    #X_N2 = 10 * x

    # plot a.1: molar fraction of each species
    plt.figure("molar fraction of species 1,2&3")
    plt.plot(x,X_O2(x),label='species1,species2')
    plt.plot(x,X_N2(x),label='species3')
    plt.legend()
    #plt.show()

    # plot a.2: molarmass of mixture
    plt.figure("molar mass of mixture")
    plt.plot(x,W_k(x))
    #plt.show()

    # plot a.3: mass fraction
    plt.figure("mass fraction of species1,2&3")
    plt.plot(x, Y_k_x1(x),label='species1')
    plt.plot(x, Y_k_x2(x), label='species1')
    plt.plot(x, Y_k_x3(x), label='species1')
    plt.legend()
    #plt.show()

    # plot a.4: density
    plt.figure("density")
    plt.plot(x, rho(x))
    #plt.show()

    # plot b: thermal conductivity
    plt.figure("mixture thermal conductivity")
    plt.plot(x, lamda(x))
    #plt.show()

    # model 1
    Model1_O2 = np.vectorize(model1_O2)(x)
    Model1_H2 = np.vectorize(model1_H2)(x)
    Model1_N2 = np.vectorize(model1_N2)(x)

    # plot c.1.a: H2 model1&2
    plt.figure("diffusion coefficients of H2 based on model 1&2")
    plt.plot(x, Model1_H2, label='model 1')
    plt.plot(x, model2_H2(x), label='model 2')
    plt.legend()
    # plot c.1.b: H2 model3&4
    plt.figure("diffusion coefficients of H2 based on model 3&4")
    plt.plot(x, model3_H2(x),label='model 3')
    plt.plot(x, model4_H2(x),label='model 4')
    plt.legend()

    # plot c.2.a: O2 model1&2
    plt.figure("diffusion coefficients of O2 based on model 1&2")
    plt.plot(x, Model1_O2, label='model 1')
    plt.plot(x, model2_O2(x), label='model 2')
    plt.legend()
    # plot c.2.b: O2 model3&4
    plt.figure("diffusion coefficients of O2 based on model 3&4")
    plt.plot(x, model3_O2(x), label='model 3')
    plt.plot(x, model4_O2(x), label='model 4')
    plt.legend()

    # plot c.3.a: N2 model1&2
    plt.figure("diffusion coefficients of N2 based on model 1&2")
    plt.plot(x, Model1_N2, label='model 1')
    plt.plot(x, model2_N2(x), label='model 2')
    plt.legend()
    # plot c.3.b: N2 model3&4
    plt.figure("diffusion coefficients of N2 based on model 3&4")
    plt.plot(x, model3_N2(x), label='model 3')
    plt.plot(x, model4_N2(x), label='model 4')
    plt.legend()
    plt.show()


part1()



