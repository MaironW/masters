# FM235 - Dinâmica de Missões Espaciais Modernas
# Author: Mairon de Souza Wolniewicz
# Date: 2025-Sep-10

# Utilities functions to be used to solve the exercises

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Generic function for Omega_x (Barcelona)
def f_Barcelona(x, mu):
    r1 = x - mu
    r2 = x - mu + 1
    sign_r1 = r1/abs(r1)
    sign_r2 = r2/abs(r2)
    return x - (1 - mu)*sign_r1/(r1**2) - mu*sign_r2/(r2**2)

# Derivative of Omega_x (Barcelona)
def df_Barcelona(x, mu):
    r1 = x - mu
    r2 = x - mu + 1
    sign_r1 = r1/abs(r1)
    sign_r2 = r2/abs(r2)
    return 1 + 2*(1 - mu)*sign_r1/(r1**3) + 2*mu*sign_r2/(r2**3)

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

# Compute the planar coordinates for the equilibrium points (x,y) in Barcelona convention
def equilibrium_points(mu, convention="Barcelona"):
    # Primaries position (Barcelona convention)
    x_P1 = mu
    x_P2 = mu - 1
    delta_x2 = -0.01 # To the left of x_P2
    delta_x3 = +0.01 # To the right of x_P2
    # Set initial guesses for the iterations based on µ and each Lagrangean point
    x1_0 = (x_P1 + x_P2)*0.5 # Between x_P1 and x_P2
    x2_0 = x_P2 + delta_x2
    x3_0 = x_P1 + delta_x3
    # For each µ, find a value for x_L1, x_L2 and x_L3
    x_L1 = newton_method(x1_0, mu, f_Barcelona, df_Barcelona)
    x_L2 = newton_method(x2_0, mu, f_Barcelona, df_Barcelona)
    x_L3 = newton_method(x3_0, mu, f_Barcelona, df_Barcelona)
    # Compute x_L4 = x_L5
    x_L4 = -0.5 + mu
    x_L5 = x_L4
    # Compute y_L4 and y_L5
    y_L4 = +3**0.5/2
    y_L5 = -y_L4
    # Define equilibrium points dict to access coordinates
    if convention=="Barcelona":
        equilibrium_points = {
            "L1" : (x_L1, 0),
            "L2" : (x_L2, 0),
            "L3" : (x_L3, 0),
            "L4" : (x_L4, y_L4),
            "L5" : (x_L5, y_L5),
            "P1" : (x_P1, 0),
            "P2" : (x_P2, 0),
        }
    else: # Caltech
        equilibrium_points = {
            "L1" : (-x_L1, 0),
            "L2" : (-x_L2, 0),
            "L3" : (-x_L3, 0),
            "L4" : (-x_L4, -y_L4),
            "L5" : (-x_L5, -y_L5),
            "P1" : (-x_P1, 0),
            "P2" : (-x_P2, 0),
        }
    return equilibrium_points

# Given x and µ, compute the Jacobi constant for the equilibrium points
def jacobi_constant(x, mu, convention="Barcelona"):
    if convention=="Barcelona":
        r1 = abs(x - mu)
        r2 = abs(x - mu + 1)
    else: # Caltech
        r1 = abs(x + mu)
        r2 = abs(x + mu - 1)
    return x**2 + 2*(1-mu)/r1 + 2*mu/r2 + mu*(1-mu)

# Return the Jacobi Constants for the entire system of equilibrium points
def jacobi_constants(equilibrium_points, mu, convention="Barcelona"):
    # Get x_L1, x_L2, x_L3
    x_L1 = equilibrium_points["L1"][0]
    x_L2 = equilibrium_points["L2"][0]
    x_L3 = equilibrium_points["L3"][0]
    # For each µ, find a value for C_L1, C_L2, C_L3
    C_L1 = jacobi_constant(x_L1, mu, convention=convention)
    C_L2 = jacobi_constant(x_L2, mu, convention=convention)
    C_L3 = jacobi_constant(x_L3, mu, convention=convention)
    # C_L4 = C_L5 is solved analytically
    C_L4 = 3
    C_L5 = 3
    # Define Jacobi Constants dict to save values
    constants = {
        "L1" : C_L1,
        "L2" : C_L2,
        "L3" : C_L3,
        "L4" : C_L4,
        "L5" : C_L5,
    }
    return constants

# Compute system µ from two body masses
def mu_from_masses(m1, m2):
    if m1<m2:
        return m1/(m1+m2)
    else:
        return m2/(m1+m2)

# Define the zero velocity limit on the planar case
def hill_curve(x, y, mu, convention="Barcelona"):
    if convention=="Barcelona":
        A = x - mu
        B = x - mu + 1
    else: # Caltech
        A = x + mu
        B = x + mu - 1
    r1 = (A**2 + y**2)**0.5
    r2 = (B**2 + y**2)**0.5
    return x**2 + y**2 + 2*(1-mu)/r1 + 2*mu/r2 + mu*(1-mu)

# Given µ, plot Hill curves for the conditions of C:
# C > C1
# C1 > C > C2
# C2 > C > C3
# C3 > C > C4
def plot_system_hill_curves(mu, convention="Barcelona"):
    # Get equilibrium points and Jacobi constants for each µ
    eq_points = equilibrium_points(mu, convention=convention)
    C_jacobi  = jacobi_constants(eq_points, mu, convention=convention)

    # Define C values to be tested
    C_cases = {
        "C>C1"    : C_jacobi["L1"]*1.01,
        "C1>C>C2" : (C_jacobi["L1"]+C_jacobi["L2"])*0.5,
        "C2>C>C3" : (C_jacobi["L2"]+C_jacobi["L3"])*0.5,
        "C3>C>C4" : (C_jacobi["L3"]+C_jacobi["L4"])*0.5,
    }

    # Create plot with equilibrium points, bodies and Hill regions
    fig, axes = plt.subplots(2,2)
    ax = axes.flatten()
    i = 0
    for C_case, C_value in C_cases.items():
        ax[i].set_title(f"µ={mu}, {C_case}")
        # Automatically plot the Lagrangian points
        for L_case, L_value in eq_points.items():
            if L_case[0] == 'L':
                ax[i].plot(L_value[0],L_value[1],'.')
                ax[i].text(L_value[0],L_value[1]-0.2,L_case, ha="center", va="center")

        # Individually plot the bodies
        ax[i].plot(eq_points["P1"][0],eq_points["P1"][1],'o', color='k', markersize=8)
        ax[i].text(eq_points["P1"][0],eq_points["P1"][1]-0.2,"P1", ha="center", va="center")
        ax[i].plot(eq_points["P2"][0],eq_points["P2"][1],'o', color='k', markersize=5)
        ax[i].text(eq_points["P2"][0],eq_points["P2"][1]-0.2,"P2", ha="center", va="center")
        ax[i].grid()
        ax[i].set_aspect("equal")

        # Fill Hill regions
        ax[i] = plot_hill_curves(ax[i], mu, C_value, convention=convention)
        i+=1

# Given a matplotlib figure axis, a µ and a Jacobi constant, plot the Hill curves
def plot_hill_curves(ax, mu, C_jacobi, convention="Barcelona"):
    grid_size = 5000

    # Define the range of the space to search for a Hill region
    x_range = np.linspace(-1.7, +1.7, grid_size)
    y_range = np.linspace(-1.7, +1.7, grid_size)
    X, Y    = np.meshgrid(x_range, y_range)

    # Fill Hill regions
    Z = hill_curve(X, Y, mu, convention=convention)
    ax.contourf(X, Y, Z, levels=[Z.min(), C_jacobi], colors=["lightgray"], alpha=0.8)
    # Draw the zero-velocity curve (boundary)
    ax.contour(X, Y, Z, levels=[C_jacobi], colors="k")
    return ax

# Build the equilibrium matrix Dxf
def linearization_matrix(x, y, z, mu, convention="Barcelona"):
    Dxf = np.zeros([6,6])

    Dxf[0,3] =  1 # df1dx4
    Dxf[1,4] =  1 # df2dx5-
    Dxf[2,5] =  1 # df3dx6
    Dxf[3,4] =  2 # df4dx5
    Dxf[4,3] = -2 # df5dx4

    # Barcelona parameters
    if convention=="Barcelona":
        A = x - mu
        B = x - mu + 1
    else: # Caltech
        A = x + mu
        B = x + mu - 1

    r1 = (A**2 + y**2 + z**2)**0.5
    r2 = (B**2 + y**2 + z**2)**0.5

    # Compute the Hessian
    Uxx = 1 - (1-mu)/r1**3 - mu/r2**3 + 3*(1-mu)*A**2/r1**5 + 3*mu*B**2/r2**5
    Uyy = 1 - (1-mu)/r1**3 - mu/r2**3 + 3*(1-mu)*y**2/r1**5 + 3*mu*y**2/r2**5
    Uzz =   - (1-mu)/r1**3 - mu/r2**3 + 3*(1-mu)*z**2/r1**5 + 3*mu*z**2/r2**5
    Uxy = 3*(1-mu)*A*y/r1**5 + 3*mu*B*y/r2**5
    Uxz = 3*(1-mu)*A*z/r1**5 + 3*mu*B*z/r2**5
    Uyz = 3*(1-mu)*y*z/r1**5 + 3*mu*y*z/r2**5

    Dxf[3:6, 0:3] = np.array([[Uxx, Uxy, Uxz],
                              [Uxy, Uyy, Uyz],
                              [Uxz, Uyz, Uzz]])
    return Dxf

# Given a list of eigenvalues and eigenvectors, return the unstable one
def get_unstable_eigenpair(values, vectors):
    reals   = np.real(values)
    pos_idx = np.argmax(reals)
    return values[pos_idx], vectors[:, pos_idx]

# General state equations for the Restricted 3 Body Problem
def state_equations(t, state, mu):
    x, y, z, vx, vy, vz = state
    A = x - mu
    B = x - mu + 1
    r1 = (A**2 + y**2 + z**2)**0.5
    r2 = (B**2 + y**2 + z**2)**0.5
    Ux = x - (1-mu) * A / r1**3 - mu * B / r2**3
    Uy = y - (1-mu) * y / r1**3 - mu * y / r2**3
    Uz = z - (1-mu) * z / r1**3 - mu * z / r2**3
    ax = +2 * vy + Ux
    ay = -2 * vx + Uy
    az = Uz
    return [vx, vy, vz, ax, ay, az]

# Find alpha algorithm as described by Parker & Chua
def find_alpha(x_eq, mu, M, eta_u, t_step=1e-2, Er = 1e-10, Ea=1e-10, alpha_min=1e-12, alpha_max=1.0):
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
        t_span = (0, M*t_step)
        solution = solve_ivp(state_equations, t_span, x_alpha, args=(mu,), rtol=1e-12, atol=1e-12)
        # Get most recent state
        px = solution.y[:, -1]
        # Linear prediction
        pLx = x_eq + alpha*np.exp(M)*eta_u
        # Check mismatch (stop condition)
        if np.linalg.norm(px - pLx) < Er*np.linalg.norm(px) + Ea:
            return alpha