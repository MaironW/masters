# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Sep-10

# 1. Cálculo dos Pontos de Equilíbrios do PR3C.

# Considere os pontos de equilíbrio colineares L1, L2, L3 do Problema Restrito de Três Corpos Circular Espacial (PR3C).
# Sendo as coordenadas destes equilíbrios representadas por (x1, 0, 0), (x2, 0, 0),
#(x3, 0, 0), respectivamente, calcule os valores de x1, x2, x3 em função do parâmetro de massa μ do modelo.
# Apresente em um único gráfico os valores de xk, k = 1, 2, 3 em função de μ para o intervalo 0 < μ < 1, juntamente com as posições dos primários P1 e P2.

import utils
import numpy as np
import matplotlib.pyplot as plt

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
    # Compute all equilibrium points for each µ
    equilibrium_points = utils.equilibrium_points(mu[i])

    # Attibute the equilibrium points x-coordinate for ease of plotting
    x_L1[i] = equilibrium_points["L1"][0]
    x_L2[i] = equilibrium_points["L2"][0]
    x_L3[i] = equilibrium_points["L3"][0]
    x_P1[i] = equilibrium_points["P1"][0]
    x_P2[i] = equilibrium_points["P2"][0]

# Plot
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
