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

import utils

# Mass of celestial bodies [kg]
body_mass = {
    "Sun"     : 1.988500e30,
    "Earth"   : 5.972190e24,
    "Jupiter" : 1.898130e27,
    "Moon"    : 7.349000e22,
    "Mars"    : 6.417100e23,
    "Venus"   : 4.867500e24,
    "Saturn"  : 5.683400e26,
    "Io"      : 8.931900e22,
    "Europa"  : 4.799800e22,
    "Titan"   : 1.345200e23,
    "Phobos"  : 1.065900e16,
    "Pluto"   : 1.303000e22,
    "Charon"  : 1.586000e21,
}

# Sort the dict by mass decrescent order
mass_sorted = dict(sorted(body_mass.items(), key=lambda item: item[1], reverse=True))
print("Solar System bodies mass:")
for b, m in mass_sorted.items():
    print(f"{b}: {m} kg")

# Compute the mass parameter µ for each system
system_mu = {
    "Earth-Moon"     : utils.mu_from_masses(body_mass["Earth"],   body_mass["Moon"]),
    "Sun-Earth"      : utils.mu_from_masses(body_mass["Sun"],     body_mass["Earth"]),
    "Sun-Venus"      : utils.mu_from_masses(body_mass["Sun"],     body_mass["Venus"]),
    "Sun-Jupiter"    : utils.mu_from_masses(body_mass["Sun"],     body_mass["Jupiter"]),
    "Sun-Saturn"     : utils.mu_from_masses(body_mass["Sun"],     body_mass["Saturn"]),
    "Jupiter-Europa" : utils.mu_from_masses(body_mass["Jupiter"], body_mass["Europa"]),
    "Saturn-Titan"   : utils.mu_from_masses(body_mass["Saturn"],  body_mass["Titan"]),
    "Jupiter-Io"     : utils.mu_from_masses(body_mass["Jupiter"], body_mass["Io"]),
    "Mars-Phobos"    : utils.mu_from_masses(body_mass["Mars"],    body_mass["Phobos"]),
    "Pluto-Charon"   : utils.mu_from_masses(body_mass["Pluto"],   body_mass["Charon"]),
}

equilibrium_points = {}
jacobi_constants   = {}
for system, mu in system_mu.items():
    # Compute all equilibrium points for each µ
    equilibrium_points[system] = utils.equilibrium_points(mu, convention="Caltech")
    x_L1 = equilibrium_points[system]["L1"][0]
    x_L2 = equilibrium_points[system]["L2"][0]
    x_L3 = equilibrium_points[system]["L3"][0]
    x_L4 = equilibrium_points[system]["L4"][0]
    x_L5 = equilibrium_points[system]["L5"][0]
    y_L4 = equilibrium_points[system]["L4"][1]
    y_L5 = equilibrium_points[system]["L5"][1]

    # Compute all Jacobi Constants for each µ
    jacobi_constants[system] = utils.jacobi_constants(equilibrium_points[system], mu, convention="Caltech")
    C_L1 = jacobi_constants[system]["L1"]
    C_L2 = jacobi_constants[system]["L2"]
    C_L3 = jacobi_constants[system]["L3"]
    C_L4 = jacobi_constants[system]["L4"]
    C_L5 = jacobi_constants[system]["L5"]

# 1. µ Table
print("\nSystem, µ")
for  system, mu in system_mu.items():
    print(f"{system},{mu}")

# 2. Coordinates Table
print("\nSystem, Coord L1, Coord L2, Coord L3, Coord L4, Coord L5")
for system, point in equilibrium_points.items():
    print(f'{system}, {point["L1"]}, {point["L2"]}, {point["L3"]}, {point["L4"]}, {point["L5"]}')

# 3. Jacobi Constants Table
print("\nSystem, C_L1, C_L2, C_L3, C_L4, C_L5")
for system, constant in jacobi_constants.items():
    print(f'{system}, {constant["L1"]}, {constant["L2"]}, {constant["L3"]}, {constant["L4"]}, {constant["L5"]}')
