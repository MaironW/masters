# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Sep-10

# 5. Estudo das Curvas de Velocidade Zero, Regiões de Hill e Regiões Inacessíveis.

# Considere o PR3C Circular Planar.
# Construa as curvas de velocidade zero para µ = 1e-4, µ = 1e-2, µ = 1e-1.
# Analise e compare as regiões acessíveis quando
# (a) C>C1,
# (b) C1>C>C2,
# (c) C2>C>C3,
# (d) C3>C>C4,C5.
# Das análises gráficas para as curvas de Hill, perceba-se que:

# (i) Para C>C2, a Região de Hill em torno de P2 é bem pequena para µ = 1e-4,
# e à medida que µ aumenta, esta região acessível aumenta.

# Ilustração de uma implicação prática de interesse real: A missão Martian Moons eXploration
# (MMX) da agência espacial japonesa JAXA em colaboração com o CNES (francês) tem lançamento
# previsto para 2026 com o objetivo de estudar as luas Phobos e Deimos de Marte.
# Planeja-se que logo após a chegada à Marte, a espaçonave será transferida para uma órbita em
# torno de Phobos.
# O semieixo maior de Phobos é da ordem de 9300 km e sua massa é bem pequena,
# de forma que o parâmetro de massa Marte-Phobos é da ordem de 1e-8.
# Sendo assim, os pontos L1 e L2 localizam-se no interior de Phobos e, portanto, uma órbita
# de estacionamento deve estar fora da chamada esfera de Hill de Phobos
# (essa esfera pode ser aproximada pela esfera de diâmetro com a distância L1 a L2).

# (ii) A região proibida às trajetórias é bem menor para µ menores que µ maiores.

# (iii) A diferença entre C1 e C2 é muito pequena para µ menores, aumentando quando µ cresce.
# Implicação prática de interesse:
# Para se realizar comparações qualitativas entre distintos sistemas do Sistema Solar,
# estes detalhes devem ser levados em conta.

import utils
import numpy as np
import matplotlib.pyplot as plt

# Define µ of interest
mu_list = [1e-4, 1e-2, 1e-1]

# Get equilibrium points and Jacobi constants for each µ
for mu in mu_list:
    equilibrium_points = utils.equilibrium_points(mu, convention="Barcelona")
    jacobi_constants = utils.jacobi_constants(equilibrium_points, mu, convention="Barcelona")
    
    plt.figure()
    plt.plot(equilibrium_points["L1"][0],equilibrium_points["L1"][1],'.')
    plt.text(equilibrium_points["L1"][0],equilibrium_points["L1"][1]-0.05,"L1", ha="center", va="center")
    plt.plot(equilibrium_points["L2"][0],equilibrium_points["L2"][1],'.')
    plt.text(equilibrium_points["L2"][0],equilibrium_points["L2"][1]-0.05,"L2", ha="center", va="center")
    plt.plot(equilibrium_points["L3"][0],equilibrium_points["L3"][1],'.')
    plt.text(equilibrium_points["L3"][0],equilibrium_points["L3"][1]-0.05,"L3", ha="center", va="center")
    plt.plot(equilibrium_points["L4"][0],equilibrium_points["L4"][1],'.')
    plt.text(equilibrium_points["L4"][0],equilibrium_points["L4"][1]+0.05,"L4", ha="center", va="center")
    plt.plot(equilibrium_points["L5"][0],equilibrium_points["L5"][1],'.')
    plt.text(equilibrium_points["L5"][0],equilibrium_points["L5"][1]-0.05,"L5", ha="center", va="center")
    plt.plot(equilibrium_points["P1"][0],equilibrium_points["P1"][1],'o', color='k', markersize=8)
    plt.text(equilibrium_points["P1"][0],equilibrium_points["P1"][1]-0.05,"P1", ha="center", va="center")
    plt.plot(equilibrium_points["P2"][0],equilibrium_points["P2"][1],'o', color='k', markersize=5)
    plt.text(equilibrium_points["P2"][0],equilibrium_points["P2"][1]-0.05,"P2", ha="center", va="center")
    plt.grid()

    # Define the regions
    C = np.linspace(jacobi_constants["L3"],jacobi_constants["L4"],100)
    for limit_C in C:
        # Solve the equation 41 for x and y

plt.show()

