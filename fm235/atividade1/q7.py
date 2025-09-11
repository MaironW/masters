# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Sep-10

# 7. Construção das variedades estável e instável dos equilíbrios colineares.

# Para o sistema Terra-Lua, a partir dos autovetores normalizados associados à direção hiperbólica,
# construa as variedades estável e instável de L1, L2 e L3.
# Para isto, sugere-se:

# (i) Apresente a aplicação do critério de escolha do parâmetro α proposto por Parker & Chua para a construção da variedade instável.
# Responda: Qual a relação entre os valores estimados de α e o módulo dos autovalores reais de cada equilíbrio colinear?

# (ii) Obtenha a variedade estável aplicando as propriedades de simetria do modelo (S1, ou S2 ou S3).

# (iii) Qual implicação numérica que você percebe existir entre a construção numérica destas variedades e as propriedades dinâmicas de cada equilíbrio?

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

# Compute the Jacobi constants for Hill regions plot
jacobi_constants = utils.jacobi_constants(equilibrium_points, mu, convention=convention)

# Compute equilibrium matrix for each case
point_list = ["L1","L2","L3","L4"]

for point in point_list:
    x = equilibrium_points[point][0]
    y = equilibrium_points[point][1]
    z = 0
    Dxf = utils.linearization_matrix(x, y, z, mu, convention="Barcelona")

    # Get eigenvalues and eigenvectors of Dxf
    eigenvalues, eigenvectors = np.linalg.eig(Dxf)
    # Sort by decreasing real part
    idx = np.argsort(-np.real(eigenvalues))
    eigenvalues  = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Get the (u)nstable eigenvalue (lambda_u) and correspondent eigenvector (eta_u)
    u_value, u_vector = utils.get_unstable_eigenpair(eigenvalues, eigenvectors)
    u_value_real = np.real(u_value) # lambda_u
    u_vector = u_vector / np.linalg.norm(u_vector) # eta_u

    # The equilibrium state is defined by [x, y, z, vx, vy, vz]
    x_eq = np.array([equilibrium_points[point][0],equilibrium_points[point][1],0,0,0,0])
    u_vector_real = np.real(u_vector)

    # Find alpha parameter if colinear equilibrium
    # And set initial condition for the integrations
    alpha = utils.find_alpha(x_eq, mu, 4, u_vector_real)
    x_alpha_plus  = x_eq + alpha * u_vector_real
    x_alpha_minus = x_eq - alpha * u_vector_real
    print(f"alpha {point}:",alpha)

    # Integration steps and duration
    t_iters = 4000
    t_range = 60

    # Forward integration (unstable)
    t_span = (0, t_range)
    t_eval = np.linspace(t_span[0], t_span[-1], t_iters)
    Wu_plus  = solve_ivp(utils.state_equations, t_span, x_alpha_plus,  t_eval=t_eval, args=(mu,), rtol=1e-10, atol=1e-12)
    Wu_minus = solve_ivp(utils.state_equations, t_span, x_alpha_minus, t_eval=t_eval, args=(mu,), rtol=1e-10, atol=1e-12)

    # Backwards integration (stable) (S1 symetry)
    t_span = (0, -t_range)
    t_eval = np.linspace(t_span[0], t_span[-1], t_iters)
    Ws_plus  = solve_ivp(utils.state_equations, t_span, x_alpha_plus,  t_eval=t_eval, args=(mu,), rtol=1e-10, atol=1e-12)
    Ws_minus = solve_ivp(utils.state_equations, t_span, x_alpha_minus, t_eval=t_eval, args=(mu,), rtol=1e-10, atol=1e-12)

    # Plot
    fig, ax = plt.subplots(1,1)
    # Orbits
    ax.plot(Wu_plus.y[0],  Wu_plus.y[1],  label="$W^{u+}$")
    ax.plot(Wu_minus.y[0], Wu_minus.y[1], label="$W^{u-}$")
    ax.plot(Ws_plus.y[0],  Ws_plus.y[1],  label="$W^{s+}$")
    ax.plot(Ws_minus.y[0], Ws_minus.y[1], label="$W^{s-}$")
    # Points
    ax.plot(equilibrium_points["P1"][0], equilibrium_points["P1"][1], "ko")
    ax.text(equilibrium_points["P1"][0], equilibrium_points["P1"][1]-0.1, "P1", ha="center", va="center")
    ax.plot(equilibrium_points["P2"][0], equilibrium_points["P2"][1], "ko")
    ax.text(equilibrium_points["P2"][0], equilibrium_points["P2"][1]-0.1, "P2", ha="center", va="center")
    for i in range(1,6):
        case = f"L{i}"
        marker = "rx" if case == point else "bx"
        ax.plot(equilibrium_points[case][0], equilibrium_points[case][1], marker)
        ax.text(equilibrium_points[case][0], equilibrium_points[case][1]-0.1, case, ha="center", va="center")
    # Hill regions
    # ax = utils.plot_hill_curves(ax, mu, jacobi_constants[point], convention=convention)

    # 2D Plot settings
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid()
    ax.legend()
    ax.axis("equal")

    # 3D Plot Position
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(Wu_plus.y[0],  Wu_plus.y[1],  Wu_plus.y[2],  label="$W^{u+}$")
    ax.plot(Wu_minus.y[0], Wu_minus.y[1], Wu_minus.y[2], label="$W^{u-}$")
    ax.plot(Ws_plus.y[0],  Ws_plus.y[1],  Ws_plus.y[2],  label="$W^{s+}$")
    ax.plot(Ws_minus.y[0], Ws_minus.y[1], Ws_minus.y[2], label="$W^{s-}$")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    # 3D Plot Velocity
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(Wu_plus.y[3],  Wu_plus.y[4],  Wu_plus.y[5],  label="$W^{u+}$")
    ax.plot(Wu_minus.y[3], Wu_minus.y[4], Wu_minus.y[5], label="$W^{u-}$")
    ax.plot(Ws_plus.y[3],  Ws_plus.y[4],  Ws_plus.y[5],  label="$W^{s+}$")
    ax.plot(Ws_minus.y[3], Ws_minus.y[4], Ws_minus.y[5], label="$W^{s-}$")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

plt.show()
