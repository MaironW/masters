# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Sep-10

# 8. Utilizando o código de sua própria autoria, se preferir,
# obtenha trajetórias do PR3C para o Sistema Terra-Lua tendo como condições iniciais as coordenadas dos 3 pontos lagrangeanos colineares e um dos triangulares.
# Apresente os gráficos destas trajetórias nos espaços xyz e ẋẏż.
# Analise e discuta os resultados obtidos, tendo em vista a análise de estabilidade linear destes pontos.

import utils
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define primary position convention
convention = "Barcelona"

# Mass of celestial bodies [kg]
body_mass = {
    "Earth"   : 5.972190e+24,
    "Moon"    : 7.349000e+22,
}

# Compute the mass parameter µ for Earth-Moon system
mu = utils.mu_from_masses(body_mass["Earth"], body_mass["Moon"])

# Compute x_eq
equilibrium_points = utils.equilibrium_points(mu, convention=convention)