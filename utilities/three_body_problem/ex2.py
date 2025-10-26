# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Oct-25

# 1. Cálculo das Órbitas Periódicas Planares do PR3C - Sistema Sol-Júpiter.

import utils
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Number of iterations
n_iter = 20

# Tolerance
tol = 1e-8

# Mass of celestial bodies [kg]
body_mass = {
    "Sun"     : 1.988500e30,
    "Jupiter" : 1.898130e27,
}

# Compute the mass parameter µ for the Sun-Jupiter system
mu = utils.mu_from_masses(body_mass["Sun"], body_mass["Jupiter"])

# Compute x_eq
equilibrium_points = utils.equilibrium_points(mu)

# Compute the Jacobi constants for Hill regions plot
jacobi_constants = utils.jacobi_constants(equilibrium_points, mu)

# Select equilibrium point
point = 'L1'

# Initial condition
dx  = -1e-2 # Small variation from position equilibrium
dvy = 0.01 # Small variation from velocity equilibrium
x0  = equilibrium_points[point][0] + dx # Initial value for x
y0  =  0.0 # Initial value for y
z0  =  0.0 # Initial value for z
vx0 =  0.0 # Initial value for vx
vy0 = 0.01 # Expected initial value for vy, to be adjusted
vz0 =  0.0 # Initial value for vz
T   =  1.1 # Expected value for the total orbit period, to be adjusted
init_state = np.array([x0,y0,z0,vx0,vy0,vz0])

# Event to detect when y crosses 0
y_event = lambda t, y, mu: y[1]
y_event.terminal  = True # Stop when event happens
y_event.direction = -1   # Search while y is decreasing

def variational_equations(t, state_M, mu):
    # Extract the state vector and the matrix M
    state  = state_M[:6]
    M_flat = state_M[6:]
    M      = M_flat.reshape((6,6))
    # Dynamic equations
    dstate = utils.state_equations(t, state, mu)
    # Jacobian matrix A(t)
    A = utils.linearization_matrix(state[0], state[1], state[2], mu,)
    # Variational equations: dM/dt = A*M
    dM = A @ M
    # Return the concatenated vector (6 + 36 variables)
    return np.concatenate((dstate, dM.flatten()))

def integrate_variational(mu, init_state, t_end):
    M0 = np.eye(6)
    y0 = np.concatenate((init_state, M0.flatten()))
    # Integrate
    t_span = (0, t_end)
    solution = solve_ivp(variational_equations, t_span, y0, args=(mu,), rtol=1e-12, atol=1e-12)
    M_flat_final = solution.y[6:,-1]
    M_final = M_flat_final.reshape((6,6))
    return M_final

def differential_correction(mu, T, state, tol=tol, n_iter=n_iter):
    for i in range(n_iter):
        # Integrate until estimated T/2 with overshooted T
        t_span = (0, T*1.5)
        solution = solve_ivp(utils.state_equations, t_span, state, args=(mu,), rtol=1e-12, atol=1e-12, events=y_event, dense_output=True)

        # Check if y=0 until estimated T/2
        if not solution.t_events[0].size:
            print(f"Convergence failed. Non-periodic orbit detected after T/2 = {T}")
            return None, T, solution

        # Find the point where y=0
        T_half = solution.t_events[0][0]
        T = 2*T_half
        solution_half = solution.sol(T_half)
        y_half  = solution_half[1]
        vx_half = solution_half[3]

        # Correction conditions
        F = np.array([y_half, vx_half])
        if np.linalg.norm(F) < tol:
            return state, T, solution # Success

        # Integrate the state transition matrix up to T/2
        M = integrate_variational(mu, state, T_half)

        # Submatrix for correction
        M_T = np.array([[M[1,0], M[1,4]],
                        [M[3,0], M[3,4]]])

        # Correction of [x0, vy0]
        delta = np.linalg.solve(M_T, -F)
        state[0] += delta[0]
        state[4] += delta[1]

    print(f"Convergence failed after {n_iter} iterations.")
    return None, T, solution

# Create Plot
fig, ax = plt.subplots(1,1)
# Points
# ax.plot(equilibrium_points["P1"][0], equilibrium_points["P1"][1], "ko")
# ax.text(equilibrium_points["P1"][0], equilibrium_points["P1"][1]-0.1, "P1", ha="center", va="center")
# ax.plot(equilibrium_points["P2"][0], equilibrium_points["P2"][1], "ko")
# ax.text(equilibrium_points["P2"][0], equilibrium_points["P2"][1]-0.1, "P2", ha="center", va="center")
# for i in range(1,6):
#     case = f"L{i}"
#     marker = "rx" if case == point else "bx"
#     ax.plot(equilibrium_points[case][0], equilibrium_points[case][1], marker)
#     ax.text(equilibrium_points[case][0], equilibrium_points[case][1]-0.1, case, ha="center", va="center")
# # Hill Regions
# ax = utils.plot_hill_curves(ax, mu, jacobi_constants[point], convention=convention)

ax.plot(equilibrium_points['L1'][0], equilibrium_points['L1'][1], 'kX')
# Build the family
state = init_state
for i in range(15):
    print(i)
    state, T, solution = differential_correction(mu, T, state)

    # Plot Orbits
    t_vals = np.linspace(0, T*0.5, 1000)
    x,y,z,vx,vy,vz = solution.sol(t_vals)

    plt.plot(x, y, label=f"Orbit {i}")
    plt.plot(x[0], y[0], 'bo')
    plt.plot(x[-1], y[-1], 'rx')

    # Stop if first divergence occured
    if state is None:
        break

    # Perturbation on velocity for next iteration
    state[4] += dvy

# 2D Plot settings
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid()
ax.legend()
ax.axis("equal")

plt.show()
