# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Sep-10

# 3. Valores de PR3C do Sistema Solar.

# Considere os seguintes subsistemas do Sistema Solar:
# Terra-Lua,
# Sol-Terra,
# Sol-Marte,
# Sol-Vênus,
# Sol-Júpiter,
# Sol-Saturno,
# Júpiter-Europa,
# Saturno-Titã,
# Júpiter-Io,
# Marte-Phobos,
# Plutão-Caronte.

# Pesquise os valores de massa disponíveis nas páginas da NASA,
# Produza uma tabela simples de duas colunas com os nomes dos corpos envolvidos nesta lista e suas massas
# Ordene essa lista de corpos em ordem decrescente de suas massas

# A partir desses valores de massas e adotando a escolha usual da Caltech para a posição dos primários,
# i.e., P1 a esquerda e P2 a direita do baricentro (origem),
# para os 11 subsistemas solicitados acima, calcule os valores de:
# (i)   parâmetro de massa µ
# (ii)  as coordenadas x dos pontos de equilíbrio L1, L2, L3, e das coordenadas x e y dos pontos de equilíbrio L4, L5
# (iii) os valores críticos da Constante de Jacobi nesses pontos de equilíbrio

# Apresente esses resultados em três tabelas:
# 1. Tabela com os nomes dos primários e os valores de µ
# 2. Tabela com 6 colunas com os nomes dos primários e os valores das coordenadas dos pontos de equilíbrio dos subsistemas solicitados
# 3. Tabela com os 4 valores críticos C1, C2, C3, e C4 = C5

# Obs.: Para obter as massas atualizadas, consulte individualmente as páginas específicas de cada corpo.
# Como por exemplo:
# https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
# https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
# que apresentam valores com mais algoritmos significativos que a tabela geral de
# https://nssdc.gsfc.nasa.gov/planetary/factsheet/

import numpy as np
import matplotlib.pyplot as plt

# Compute the system µ based on each body mass
def system_mu(m1, m2):
    if m1>m2:
        return m1/(m1 + m2)
    else:
        return m2/(m1 + m2)

# Generic function for Omega_x (Caltech)
def f_Caltech(x, mu):
    A = x + mu
    B = x + mu - 1
    sign_A = A/abs(A)
    sign_B = B/abs(B)
    return x - (1 - mu)*sign_A/(A**2) - mu*sign_B/(B**2)

# Derivative of Omega_x (Caltech)
def df_Caltech(x, mu):
    A = x + mu
    B = x + mu - 1
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
def jacobi_constants(x, mu):
    A = abs(x + mu)
    B = abs(x + mu - 1)
    return x**2 + 2*(1-mu)/A + 2*mu/B + mu*(1-mu)

# Compute the planar coordinates for the equilibrium points (x,y)
def equilibrium_coordinates(mu):
    # Primaries position (Caltech convention)
    x_P1 = mu
    x_P2 = mu - 1
    delta_x2 = +0.01 # To the right of x_P2
    delta_x3 = -0.01 # To the left of x_P2

    # Set initial guesses for the iterations based on µ and each Lagrangean point
    x1_0 = (x_P1 + x_P2)*0.5 # Between x_P1 and x_P2
    x2_0 = x_P2 + delta_x2
    x3_0 = x_P1 + delta_x3

    # For each µ, find a value for x_L1, x_L2 and x_L3
    x_L1 = newton_method(x1_0, mu, f_Caltech, df_Caltech)
    x_L2 = newton_method(x2_0, mu, f_Caltech, df_Caltech)
    x_L3 = newton_method(x3_0, mu, f_Caltech, df_Caltech)

    # Compute x_L4 = x_L5
    x_L4 = -0.5 + mu
    x_L5 = x_L4

    # Compute y_L4 and y_L5
    y_L4 = +3**0.5/2
    y_L5 = -x_L4

    equilibrium_points = {
        "L1" : (x_L1, 0),
        "L2" : (x_L2, 0),
        "L3" : (x_L3, 0),
        "L4" : (x_L4, y_L4),
        "L5" : (x_L5, y_L5),
    }

    return equilibrium_points


# Mass of celestial bodies [kg]
body_mass = {
    "Sun"     : 1.988500e+30,
    "Earth"   : 5.972190e+24,
    "Jupiter" : 1.898130e+27,
    "Moon"    : 7.349000e+22,
    "Mars"    : 6.417100e+23,
    "Venus"   : 4.867500e+24,
    "Saturn"  : 5.683400e+26,
    "Io"      : 8.931900e+22,
    "Europa"  : 4.799800e+22,
    "Titan"   : 1.345200e+23,
    "Phobos"  : 1.065900e+16,
    "Pluto"   : 1.303000e+22,
    "Charon"  : 1.586000e+21,
}

# Sort the dict by mass decrescent order
mass_sorted = dict(sorted(body_mass.items(), key=lambda item: item[1], reverse=True))
print("Solar System bodies mass:")
for b, m in mass_sorted.items():
    print(f"{b}: {m} kg")

# Compute the mass parameter µ for each system
system_mu = {
    "Earth-Moon"     : system_mu(body_mass["Earth"],   body_mass["Moon"]),
    "Sun-Earth"      : system_mu(body_mass["Sun"],     body_mass["Earth"]),
    "Sun-Venus"      : system_mu(body_mass["Sun"],     body_mass["Venus"]),
    "Sun-Jupiter"    : system_mu(body_mass["Sun"],     body_mass["Jupiter"]),
    "Sun-Saturn"     : system_mu(body_mass["Sun"],     body_mass["Saturn"]),
    "Jupiter-Europa" : system_mu(body_mass["Jupiter"], body_mass["Europa"]),
    "Saturn-Titan"   : system_mu(body_mass["Saturn"],  body_mass["Titan"]),
    "Jupiter-Io"     : system_mu(body_mass["Jupiter"], body_mass["Io"]),
    "Mars-Phobos"    : system_mu(body_mass["Mars"],    body_mass["Phobos"]),
    "Pluto-Charon"   : system_mu(body_mass["Pluto"],   body_mass["Charon"]),
}

equilibrium_points = {}
jacobi_constants   = {}
for system, mu in system_mu.items():
    equilibrium_points[system] = equilibrium_coordinates(mu)