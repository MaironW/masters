# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz

# 2. Cálculo dos Valores Críticos da Constante de Jacobi.

# Calcule o valor da constante de Jacobi de cada um dos cinco pontos de equilíbrio do PR3C em função do parâmetro de massa.
# Apresente em um único gráfico os valores de C1 , C2 , C3 , C4 = C5 como função de µ, para µ ∈ (0, 1).
# Definindo Ck = −2 Ek , plote Ek × µ, para k = 1, 2, 3, 4, 5.

import numpy as np
import matplotlib.pyplot as plt

# Generic function for Omega_x
def f(x, mu):
    A = x - mu
    B = x - mu + 1
    sign_A = A/abs(A)
    sign_B = B/abs(B)
    return x - (1 - mu)*sign_A/(A**2) - mu*sign_B/(B**2)

# Derivative of Omega_x
def df(x, mu):
    A = x - mu
    B = x - mu + 1
    sign_A = A/abs(A)
    sign_B = B/abs(B)
    return 1 + 2*(1 - mu)*sign_A/(A**3) + 2*mu*sign_B/(B**3)

# Newton-Raphson method to solve non-linear equations
def newton_method(x_0, mu, f, df, tol=1e-12, max_iter=10000):
    i     = 0
    x_n   = x_0
    error = 1000
    while (error > tol) and (i < max_iter):
        x_n1 = x_n - f(x_n, mu)/df(x_n, mu)
        error = abs(x_n1 - x_n)
        x_n = x_n1
        i += 1
    return x_n1

# Given x and µ, compute the Jacobi constant for the equilibrium points
def jacobi_constant(x, mu):
    A = abs(x-mu)
    B = abs(x-mu+1)
    return x**2 + 2*(1-mu)/A + 2*mu/B + mu*(1-mu)

# Number of iterations
N = 100

# Variables to save data
x_L1 = np.zeros(N)
x_L2 = np.zeros(N)
x_L3 = np.zeros(N)
x_P1 = np.zeros(N)
x_P2 = np.zeros(N)
C_L1 = np.zeros(N)
C_L2 = np.zeros(N)
C_L3 = np.zeros(N)
C_L4 = np.ones(N)*3 # Solved analytically
C_L5 = np.ones(N)*3 # Solved analytically

# Iterate µ from 0 to 1 (with tolerances to avoid singularities)
mu = np.linspace(0+1e-12, 1-1e-12, N)

for i in range(N):
    # Primaries position (following Barcelona convention)
    x_P1[i] = mu[i]
    x_P2[i] = -1 + mu[i]

    # Set initial guesses for the iterations based on µ and each Lagrangean point
    x1_0 = (x_P1[i] + x_P2[i])*0.5 # Between x_P1 and x_P2
    x2_0 = x_P2[i] - 0.01   # To the left of x_P2
    x3_0 = x_P1[i] + 0.01   # To the right of x_P2

    # For each µ, find a value for x_L1, x_L2 and x_L3
    x_L1[i] = newton_method(x1_0, mu[i], f, df)
    x_L2[i] = newton_method(x2_0, mu[i], f, df)
    x_L3[i] = newton_method(x3_0, mu[i], f, df)

    # For each µ, find a value for C_L1, C_L2, C_L3
    C_L1[i] = jacobi_constant(x_L1[i], mu[i])
    C_L2[i] = jacobi_constant(x_L2[i], mu[i])
    C_L3[i] = jacobi_constant(x_L3[i], mu[i])

plt.figure("CR3BP Jacobi Constants (Barcelona et al convention)")
plt.grid()
plt.plot(mu, -C_L1/2.0, label="$E_{L1}$")
plt.plot(mu, -C_L2/2.0, label="$E_{L2}$")
plt.plot(mu, -C_L3/2.0, label="$E_{L3}$")
plt.plot(mu, -C_L4/2.0, label="$E_{L4}=E_{L5}$")
plt.xlabel("$\mu$")
plt.ylabel("$C_k$")
plt.legend()
plt.show()
