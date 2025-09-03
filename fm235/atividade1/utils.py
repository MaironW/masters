import numpy as np
import matplotlib.pyplot as plt

# Generic function for Omega_x (Barcelona)
def f_Barcelona(x, mu):
    A = x - mu
    B = x - mu + 1
    sign_A = A/abs(A)
    sign_B = B/abs(B)
    return x - (1 - mu)*sign_A/(A**2) - mu*sign_B/(B**2)

# Derivative of Omega_x (Barcelona)
def df_Barcelona(x, mu):
    A = x - mu
    B = x - mu + 1
    sign_A = A/abs(A)
    sign_B = B/abs(B)
    return 1 + 2*(1 - mu)*sign_A/(A**3) + 2*mu*sign_B/(B**3)

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
        A = abs(x - mu)
        B = abs(x - mu + 1)
    else: # Caltech
        A = abs(x + mu)
        B = abs(x + mu - 1)
    return x**2 + 2*(1-mu)/A + 2*mu/B + mu*(1-mu)

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