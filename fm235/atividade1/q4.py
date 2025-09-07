# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Sep-10

# 4. O Problema de Copenhagen.

# O caso especial do PR3C com µ = 0.5 é conhecido como o Problema de Copenhagen e representa o caso simétrico em que as massas dos primários são iguais. 
# Perceba e comente dois aspectos relativos a esse caso:
# (i) No caso particular de µ = 0.5, há uma simetria no espaço de fases envolvendo a posição dos pontos lagrangeanos e respectivos valores de constante de Jacobi (conforme gráficos solicitados nas Questões 1 e 2). 
# Comente as possibilidades de transporte em função de C no Problema de Copenhagen.
# (ii) O Problema de Copenhagen pode ser considerado o caso de limite superior para o valor de µ, devido a outra simetria que existe em torno de µ = 0.5 no plano paramétrico (correspondendo a trocar P1 por P2).

import utils

# Define primary position convention
convention = "Barcelona"

# Define µ as 0.5
mu = 0.5

# Compute the coordinates of the equilibrium points
equilibrium_points = utils.equilibrium_points(mu, convention=convention)

# Output equilibrium points
print("Equilibrium Points Coordinates:")
print(equilibrium_points)

# Compute all Jacobi Constants
jacobi_constants = utils.jacobi_constants(equilibrium_points, mu, convention=convention)

# Output Jacobi constants
print("Jacobi Constants:")
print(jacobi_constants)



