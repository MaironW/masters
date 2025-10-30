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
def integrate_variational(mu, init_state, t_end, terminate=True):
    def y_event(t, y, mu):
        return y[1]
    y_event.terminal = terminate
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
        solution = integrate_variational(mu, state, T_half, terminate=True)
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
                break

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
    "Jupiter" : 1.898130e27,
}

# Compute the mass parameter µ for the Sun-Jupiter system
mu = utils.mu_from_masses(body_mass["Sun"], body_mass["Jupiter"])

# Compute equilibrium points
equilibrium_points = utils.equilibrium_points(mu)

# Get solutions
def compute_solutions():
    # Lists to store solutions
    init_state_list = []
    tf = 1.5
    for x0 in np.linspace(equilibrium_points['L1'][0]-0.002, -0.99, 20):
        # Brute force first guesses for vy0
        vy0_new = None
        for vy0 in np.linspace(0.001, 0.6, 20):
            # vy0_new will be None until a valid solution is found.
            # In this case, keep using the brute force guesses.
            if vy0_new is not None:
                vy0 = vy0_new
            init_state = np.array([x0, 0, 0, 0, vy0, 0])
            solution, T_half = differential_correction(mu, init_state, tf)

            if solution.t_events[0].size > 0:
                T_half = solution.t_events[0][0]
                t_vals = np.linspace(0, T_half*2, 1000)
                state = solution.sol(t_vals)[:6]
                # If xf > x0, save updated init state
                if state[0][-1] > state[0][0]:
                    x,y,z,vx,vy,vz = solution.sol(0)[:6]
                    print(x,y,z,vx,vy,vz,T_half)
                    vy0_new = vy
                    init_state_list.append([x,y,z,vx,vy,vz,T_half])
                    tf = T_half*2
                    break # Stop iterating the velocity if solution found
    return init_state_list

# Load solutions from file in format [x,y,z,vx,vy,vz,T_half]
def load_solutions(file):
    init_state_list = np.loadtxt(file, delimiter=",")
    return init_state_list

# Integrate a full orbit given initial states
def integrate_solutions(init_state_list):
    full_state_list = []
    for init_state in init_state_list:
        x,y,z,vx,vy,vz,T_half = init_state
        tf = T_half
        solution = integrate_variational(mu, init_state[:6], 2*tf, terminate=False)
        t_vals = np.linspace(0, 2*tf, 1000)
        full_state = solution.sol(t_vals)[:6]
        full_state_list.append(full_state)
    return full_state_list

# init_state_list = compute_solutions()
# np.savetxt("lyapunov.txt", init_state_list, delimiter=",")

init_state_list = load_solutions("lyapunov.txt")
full_state_list = integrate_solutions(init_state_list)

# Process results
x0_list     = []
vy0_list    = []
T_half_list = []
jacobi_list = []
for state in init_state_list:
    x0     = state[0]
    vy0    = state[4]
    T_half = state[6]
    x0_list.append(x0)
    vy0_list.append(vy0)
    T_half_list.append(T_half)
    jacobi_list.append(utils.jacobi_constant(x0, mu))

# i) Orbit plot
fig, ax = plt.subplots(1,2)
ax[0].plot(equilibrium_points['L1'][0], equilibrium_points['L1'][1], 'kX')
ax[0].plot(equilibrium_points['L2'][0], equilibrium_points['L2'][1], 'kX')
for state in full_state_list:
    ax[0].plot(state[0],state[1], label=f"x0: {state[0][0]:.3g}, vy0: {state[4][0]:.3g}, tf: {T_half:.3g}")
    ax[1].plot(state[3],state[4], label=f"x0: {state[0][0]:.3g}, vy0: {state[4][0]:.3g}, tf: {T_half:.3g}")
ax[0].set_xlabel("x")
ax[0].set_ylabel("y")
ax[0].grid()
ax[0].legend()
ax[0].axis("equal")
ax[1].set_xlabel("x_dot")
ax[1].set_ylabel("y_dot")
ax[1].grid()
ax[1].legend()
ax[1].axis("equal")

# ii)
fig, ax = plt.subplots(2,1)
ax[0].plot(x0_list, vy0_list)
ax[1].plot(x0_list, T_half_list)
ax[1].set_xlabel("x_0")
ax[0].set_ylabel("y_dot_0")
ax[1].set_ylabel("T_half")

# C)
fig, ax = plt.subplots(2,1)
ax[0].plot(x0_list, jacobi_list)
ax[1].plot(x0_list, np.array(T_half_list)*2)
ax[1].set_xlabel("x_0")
ax[0].set_ylabel("C")
ax[1].set_ylabel("T")

plt.show()
