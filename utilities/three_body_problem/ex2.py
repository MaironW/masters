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
dx  = 5e-4 # Small variation from position equilibrium
x0  = equilibrium_points[point][0] + dx # Initial value for x
y0  =  0.0 # Initial value for y
z0  =  0.0 # Initial value for z
vx0 =  0.0 # Initial value for vx
vy0 =  0.05 # Expected initial value for vy, to be adjusted
vz0 =  0.0 # Initial value for vz
T_half = 1.1 # Expected value for the half orbit period, to be adjusted
state  = np.array([x0,y0,z0,vx0,vy0,vz0]) # Initial state

# Detect when there is a crossing
def find_y_crossing(solution):
    y_vals = solution.y[1]
    t_vals = solution.t
    for i in range(len(y_vals)-1):
        if y_vals[i]*y_vals[i+1] < 0:
            t_cross = t_vals[i] - y_vals[i]*(t_vals[i+1]-t_vals[i])/(y_vals[i+1]-y_vals[i])
            return t_cross
    return None

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

# Integrate the state transition matrix up to T/2
def integrate_variational(mu, init_state, t_end):
    # M =  dx/dx  dx/dy  dx/dz  dx/dvx  dx/dvy    dx/dvz
    #      dy/dx  dy/dy  dy/dz  dy/dvx [dy/dvy]   dy/dvz
    #      dz/dx  dz/dy  dz/dz  dz/dvx  dz/dvy    dz/dvz
    #     dvx/dx dvx/dy dvx/dz dvx/dvx [dvx/dvy] dvx/dvz
    #     dvy/dx dvy/dy dvy/dz dvy/dvx  dvy/dvy  dvy/dvz
    #     dvz/dx dvz/dy dvz/dz dvz/dvx  dvz/dvy  dvz/dvz
    M0 = np.eye(6).flatten()
    y0 = np.concatenate((init_state, M0))
    t_span = (0, t_end)
    solution = solve_ivp(variational_equations, t_span, y0, args=(mu,), rtol=1e-12, atol=1e-12, dense_output=True)
    state_final  = solution.y[:6,-1]
    M_flat_final = solution.y[6:,-1]
    M_final = M_flat_final.reshape((6,6))
    return state_final, M_final, solution

def differential_correction(mu, T_half, state, tol=tol, n_iter=n_iter):
    for i in range(n_iter):
        # Integrate until estimated T/2 with overshooted T
        t_span = (0, T_half*2)
        # solution = solve_ivp(utils.state_equations, t_span, state, args=(mu,), rtol=1e-12, atol=1e-12, dense_output=True)
        state, M, solution = integrate_variational(mu, state, T_half*2)

        # # Find the point where y=0
        # T_crossing = find_y_crossing(solution)
        # if T_crossing is None:
        #     print(f"Convergence failed. No y=0 crossing after T/2 = {T_half}")
        #     state_half = solution.sol(T_half)
        #     return state_half, T_half, solution
        # T_half = T_crossing

        # Get useful variables at half period
        state_half = state
        deriv_half = np.array(utils.state_equations(T_half, state_half, mu))
        y_half  = state_half[1]
        vx_half = state_half[3]
        vy_half = state_half[4]
        ax_half = deriv_half[3]

        # Correction conditions
        F = np.array([y_half, vx_half])
        if np.linalg.norm(F) < tol:
            return state, T_half, solution # Success

        # Get useful parameters from the state transition matrix at half period
        dy_dvy  = M[1,4]
        dvx_dvy = M[3,4]

        # Build the minimal state transition matrix
        M_sub = np.array([[ dy_dvy, vy_half],
                          [dvx_dvy, ax_half]])

        # Correction of [vy0, tf]
        delta = np.linalg.solve(M_sub, -F)
        state[4] += delta[0]
        T_half   += delta[1]

    print(f"Convergence failed after {n_iter} iterations.")
    return state, T_half, solution

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
for i in range(10):
    print(i)
    state, T_half, solution = differential_correction(mu, T_half, state)

    # Plot Orbits
    t_vals = np.linspace(0, T_half, 1000)
    print(solution.sol(t_vals))
    x,y,z,vx,vy,vz = solution.sol(t_vals)[:6]

    plt.plot(x, y, label=f"Orbit {i}")
    plt.plot(x[0], y[0], 'bo')
    plt.plot(x[-1], y[-1], 'rx')

    # Stop if first divergence occured
    # if state is None:
        # break

    # Increment to get family
    state[0] += dx

# 2D Plot settings
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid()
ax.legend()
ax.axis("equal")

fig, ax = plt.subplots(1,1)
ax.plot(t_vals, x)
ax.plot(t_vals, y)
ax.plot(t_vals, z)

fig, ax = plt.subplots(1,1)
ax.plot(t_vals, vx)
ax.plot(t_vals, vy)
ax.plot(t_vals, vz)

plt.show()
