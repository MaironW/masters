import utils
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def variational_equations(t, state_M, mu):
    # Extract the state vector and the matrix M
    state  = state_M[:6]
    M_flat = state_M[6:]
    M      = M_flat.reshape((6,6))
    # Dynamic equations
    dstate = utils.state_equations(t, state, mu)
    # Jacobian matrix A(t)
    A = utils.linearization_matrix(state[0], state[1], state[2], mu)
    # Variational equations: dM/dt = A*M
    dM = A @ M
    # Return the concatenated vector (6 + 36 variables)
    return np.concatenate((dstate, dM.flatten()))

# Integrate the state transition matrix up to T/2, when the solution crosses y=0
def integrate_variational(mu, init_state, t_end):
    def y_event(t, y, mu):
        return y[1]
    y_event.terminal = True
    y_event.direction = -1

    # M =  dx/dx  dx/dy  dx/dz  dx/dvx  dx/dvy    dx/dvz
    #      dy/dx  dy/dy  dy/dz  dy/dvx [dy/dvy]   dy/dvz
    #      dz/dx  dz/dy  dz/dz  dz/dvx  dz/dvy    dz/dvz
    #     dvx/dx dvx/dy dvx/dz dvx/dvx [dvx/dvy] dvx/dvz
    #     dvy/dx dvy/dy dvy/dz dvy/dvx  dvy/dvy  dvy/dvz
    #     dvz/dx dvz/dy dvz/dz dvz/dvx  dvz/dvy  dvz/dvz
    M0 = np.eye(6).flatten()
    y0 = np.concatenate((init_state, M0))
    t_span = (0, t_end)
    solution = solve_ivp(variational_equations, t_span, y0, args=(mu,), rtol=1e-12, atol=1e-12, events=y_event, dense_output=True)
    return solution

def differential_correction(mu, init_state, t_end, tol=1e-8, n_iter=10):
    # Make a copy of the initial state to not alter it
    state = np.copy(init_state)
    T_half = t_end
    for i in range(n_iter):
        # Integrate the augmented state
        solution = integrate_variational(mu, state, T_half)
        # The solution may be valid if the solution crosses the y=0 axis
        if solution.t_events[0].size > 0:
            t_end = solution.t_events[0][0]
            state_half = solution.sol(t_end)[:6]
            deriv_half = utils.state_equations(t_end, state_half, mu)
            y_half  = state_half[1]
            vx_half = state_half[3]
            vy_half = state_half[4]
            ax_half = deriv_half[3]

            # Correction condition
            F = np.array([y_half, vx_half])
            if np.linalg.norm(F) < tol:
                return solution, T_half

            # Get useful parameters from the state transition matrix at half period
            M = solution.sol(t_end)[6:].reshape((6,6))
            dy_dvy  = M[1,4]
            dvx_dvy = M[3,4]

            # Build the minimal state transition matrix
            M_sub = np.array([[ dy_dvy, vy_half],
                              [dvx_dvy, ax_half]])

            # Correction of [vy0, tf]
            delta = np.linalg.solve(M_sub, -F)
            state[4] += delta[0]
            T_half   += delta[1]

    return solution, T_half


# Mass of celestial bodies [kg]
body_mass = {
    "Sun"     : 1.988500e30,
    "Earth"   : 5.972190e24,
    "Jupiter" : 1.898130e27,
    "Moon"    : 7.349000e22,
}

# Compute the mass parameter µ for the Earth-Moon system
mu = utils.mu_from_masses(body_mass["Earth"], body_mass["Moon"])

# Compute equilibrium points
equilibrium_points = utils.equilibrium_points(mu)

# Compute the Jacobi constants for Hill regions plot
jacobi_constants = utils.jacobi_constants(equilibrium_points, mu)

# Create Plot
fig, ax = plt.subplots(1,1)
ax.plot(equilibrium_points['L1'][0], equilibrium_points['L1'][1], 'kX')
ax.plot(equilibrium_points['L2'][0], equilibrium_points['L2'][1], 'kX')

# Get solutions
tf = 1.5
for dx in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]:
    x0 = equilibrium_points['L1'][0] - dx
    print(x0)
    # Brute force first guesses for vy0
    vy0_new = None
    for vy0 in np.linspace(0.07, dx*10, 10):
        init_state = np.array([x0, 0, 0, 0, vy0, 0])
        # solution = integrate_variational(mu, init_state, tf)
        solution, T_half = differential_correction(mu, init_state, tf)

        if solution.t_events[0].size > 0:
            T_half = solution.t_events[0][0]
            t_vals = np.linspace(0, T_half*2, 1000)
            state = solution.sol(t_vals)[:6]
            # If xf > x0, save updated init state
            if state[0][-1] > state[0][0]:
                init_state = solution.sol(0)
                solution = solve_ivp(variational_equations, (0,2*tf), init_state, args=(mu,), rtol=1e-12, atol=1e-12, dense_output=True)
                full_state = solution.sol(t_vals)[:6]
                ax.plot(full_state[0],full_state[1], label=f"vy0: {init_state[4]}, tf: {T_half}")
                break # Stop iterating the velocity if solution found

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid()
ax.legend()
ax.axis("equal")
plt.show()
