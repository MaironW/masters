# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz

# 1. Cálculo dos Pontos de Equilíbrios do PR3C.

# Considere os pontos de equilíbrio colineares L1, L2, L3 do Problema Restrito de Três Corpos Circular Espacial (PR3C).
# Sendo as coordenadas destes equilíbrios representadas por (x1, 0, 0), (x2, 0, 0),
#(x3, 0, 0), respectivamente, calcule os valores de x1, x2, x3 em função do parâmetro de massa μ do modelo.
# Apresente em um único gráfico os valores de xk, k = 1, 2, 3 em função de μ para o intervalo 0 < μ < 1, juntamente com as posições dos primários P1 e P2.

import numpy as np
import matplotlib.pyplot as plt

# Generic function for Omega_x
def f(x, mu):
    A = x - mu
    B = x + 1 - mu
    sign_A = A/abs(A)
    sign_B = B/abs(B)
    return x - (1 - mu)*sign_A/(A**2) - mu*sign_B/(B**2)

# Derivative of Omega_x
def df(x, mu):
    A = x - mu
    B = x + 1 - mu
    sign_A = A/abs(A)
    sign_B = B/abs(B)
    return 1 + 2*(1 - mu)*sign_A/(A**3) + 2*sign_B/(B**3)

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

# Number of iterations
N = 100

# Variables to save data
x_L1 = np.zeros(N)
x_L2 = np.zeros(N)
x_L3 = np.zeros(N)
x_P1 = np.zeros(N)
x_P2 = np.zeros(N)

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

plt.figure("CR3BP Equilibrium points (Barcelona et al convention)")
plt.grid()
plt.plot(mu, x_L1, label="$x_{L1} (x_{P2} < x_{L1} < x_{P1})$")
plt.plot(mu, x_L2, label="$x_{L2} (x_{L2} < x_{P2})$")
plt.plot(mu, x_L3, label="$x_{L3} (x_{P1} < x_{L3})$")
plt.plot(mu, x_P1, label="$x_{P1}$")
plt.plot(mu, x_P2, label="$x_{P2}$")
plt.xlabel("$\mu$")
plt.ylabel("$x_k$")
plt.legend()
plt.show()
