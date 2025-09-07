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

def state_equations(t, state):
    x, y, vx, vy = state
    A = x - mu
    B = x - mu + 1
    r1 = (A**2 + y**2)**0.5
    r2 = (B**2 + y**2)**0.5
    Ux = x - (1-mu) * A / r1**3 - mu * B / r2**3
    Uy = y - (1-mu) * y / r1**3 - mu * y / r2**3
    ax = +2 * vy + Ux
    ay = -2 * vx + Uy
    return [vx, vy, ax, ay]

# Mass of celestial bodies [kg]
body_mass = {
    "Earth"   : 5.972190e+24,
    "Moon"    : 7.349000e+22,
}

# Compute the mass parameter µ for Earth-Moon system
mu = utils.mu_from_masses(body_mass["Earth"], body_mass["Moon"])

# Compute x_eq
equilibrium_points = utils.equilibrium_points(mu, convention=convention)

# Create a dict to store the eigenvalues and eigenvectors for each case
eigen_data = {}

# Compute equilibrium matrix for each case
point_list = ["L1","L2","L3"]

for point in point_list:
    x = equilibrium_points[point][0]
    y = equilibrium_points[point][1]
    z = 0
    Dxf = utils.linearization_matrix(x, y, z, mu, convention="Barcelona", planar=True)

    # Get eigenvalues and eigenvectors of Dxf
    eigenvalues, eigenvectors = np.linalg.eig(Dxf)

    # Get the (u)nstable eigenvalue (lambda_u) and correspondent eigenvector (eta_u)
    u_value, u_vector = utils.get_unstable_eigenpair(eigenvalues, eigenvectors)
    u_value_real = np.real(u_value) # lambda_u
    u_vector = u_vector / np.linalg.norm(u_vector) # eta_u

    # The equilibrium state is defined by [x, y, vx, vy]
    x_eq = np.array([equilibrium_points[point][0],0,0,0])
    u_vector_real = np.real(u_vector)

    # Find alpha algorithm as described by Parker & Chua
    def find_alpha(x_eq, M, eta_u, propagate, Er = 1e-12, Ea=1e-12, alpha_min=1e-12, alpha_max=1.0):
        alpha = 2*alpha_max
        while 1:
            # Halve alpha every iteration
            alpha /= 2.0
            # Lower limit alpha by alpha_min (last stop condition)
            if alpha <= alpha_min:
                return alpha_min
            # Perturbed initial condition
            x_alpha = x_eq + alpha*eta_u
            # Nonlinear propagation
            px = propagate(x_alpha, M)
            # Linear prediction
            pLx = x_eq + alpha*np.exp(M)*eta_u
            # Check mismatch (stop condition)
            if np.linalg.norm(px - pLx) < Er*np.linalg.norm(px) + Ea:
                return alpha

    # Propagation function to find alpha given M steps
    def propagate(x0, M, t_step=1e-2):
        t_span = (0, M*t_step)
        solution = solve_ivp(state_equations, t_span, x0, rtol=1e-12, atol=1e-12)
        # Return the last state
        return solution.y[:, -1]

    # Find alpha parameter
    alpha = find_alpha(x_eq, 4, u_vector_real, propagate=propagate)

    # Set initial condition for the integrations
    x_alpha_plus  = x_eq + alpha * u_vector_real
    x_alpha_minus = x_eq - alpha * u_vector_real

    # Integration steps and duration
    t_iters = 10000
    t_range = 60

    # Forward integration (unstable)
    t_span = (0, t_range/abs(u_value_real))
    t_eval = np.linspace(t_span[0], t_span[-1], t_iters)
    sol_plus_f  = solve_ivp(state_equations, t_span, x_alpha_plus,  t_eval=t_eval, rtol=1e-10, atol=1e-12)
    sol_minus_f = solve_ivp(state_equations, t_span, x_alpha_minus, t_eval=t_eval, rtol=1e-10, atol=1e-12)

    # Backwards integration (stable)
    t_span = (0, -t_range/abs(u_value_real))
    t_eval = np.linspace(t_span[0], t_span[-1], t_iters)
    sol_plus_b  = solve_ivp(state_equations, t_span, x_alpha_plus,  t_eval=t_eval, rtol=1e-10, atol=1e-12)
    sol_minus_b = solve_ivp(state_equations, t_span, x_alpha_minus, t_eval=t_eval, rtol=1e-10, atol=1e-12)

    plt.figure()
    plt.plot(sol_plus_f.y[0],  sol_plus_f.y[1],  label="$W^{u+}$")
    plt.plot(sol_minus_f.y[0], sol_minus_f.y[1], label="$W^{u-}$")
    plt.plot(sol_plus_b.y[0],  sol_plus_b.y[1],  label="$W^{s+}$")
    plt.plot(sol_minus_b.y[0], sol_minus_b.y[1], label="$W^{s-}$")
    plt.plot(equilibrium_points["P1"][0], equilibrium_points["P1"][1], "ko")
    plt.text(equilibrium_points["P1"][0], equilibrium_points["P1"][1]-0.1, "P1", ha="center", va="center")
    plt.plot(equilibrium_points["P2"][0], equilibrium_points["P2"][1], "ko")
    plt.text(equilibrium_points["P2"][0], equilibrium_points["P2"][1]-0.1, "P2", ha="center", va="center")
    for i in range(1,6):
        case = f"L{i}"
        marker = "rx" if case == point else "bx"
        plt.plot(equilibrium_points[case][0], equilibrium_points[case][1], marker)
        plt.text(equilibrium_points[case][0], equilibrium_points[case][1]-0.1, case, ha="center", va="center")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.legend()
    plt.axis("equal")
    plt.tight_layout()

    plt.show()
