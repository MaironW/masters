# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Sep-10

# 2. Cálculo dos Valores Críticos da Constante de Jacobi.

# Calcule o valor da constante de Jacobi de cada um dos cinco pontos de equilíbrio do PR3C em função do parâmetro de massa.
# Apresente em um único gráfico os valores de C1 , C2 , C3 , C4 = C5 como função de µ, para µ ∈ (0, 1).
# Definindo Ck = −2 Ek , plote Ek × µ, para k = 1, 2, 3, 4, 5.

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
C_L1 = np.zeros(N)
C_L2 = np.zeros(N)
C_L3 = np.zeros(N)
C_L4 = np.ones(N)*3 # Solved analytically
C_L5 = np.ones(N)*3 # Solved analytically

# Iterate µ from 0 to 1 (with tolerances to avoid singularities)
mu = np.linspace(0+1e-12, 1-1e-12, N)

for i in range(N):
    # Compute all equilibrium points for each µ
    equilibrium_points = utils.equilibrium_points(mu[i], convention="Barcelona")

    # Attibute the equilibrium points x-coordinate for ease of plotting
    x_L1[i] = equilibrium_points["L1"][0]
    x_L2[i] = equilibrium_points["L2"][0]
    x_L3[i] = equilibrium_points["L3"][0]

    # For each µ, find a value for C_L1, C_L2, C_L3
    C_L1[i] = utils.jacobi_constant(x_L1[i], mu[i], convention="Barcelona")
    C_L2[i] = utils.jacobi_constant(x_L2[i], mu[i], convention="Barcelona")
    C_L3[i] = utils.jacobi_constant(x_L3[i], mu[i], convention="Barcelona")

# Plot
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
