import numpy as np

# Parameters for Module NAV_UKF

# Process noise mapping
G = np.zeros((6, 3))
G[3:6, :] = np.eye(3)

# Acceleration noise level
sigma_a = 1e-6

# Process noise covariance
Q = sigma_a**2 * np.eye(3)

NAV_UKF_par = {
    "dt"        : 600, # [s] Integration step. Update according to the SIM_dt
    "n_states"  : 6,   # Number of estimated states [pos[3], vel[3]]
    "G"         : G,   # Process noise mapping
    "Q"         : Q,   # Process noise covariance
    "alpha"     : 1e-3,
    "beta"      : 2,
    "kappa"     : 0,

    # Initial values for the UKF
    "x_est_ini" : np.array([1.62228612e+08, 2.42235379e+07, 1.78903662e+07, -0.32883677, 27.86282899, 13.43536119]) + np.array([10, 10, 10, 0.01, 0.01, 0.01]), # Initial state
    "P_ini"     : np.ones((6,6)), # Initial covariance
    "y_ini"     : 0, # Initial innovation
}
