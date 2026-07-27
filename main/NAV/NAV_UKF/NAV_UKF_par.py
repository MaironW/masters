import numpy as np

# Parameters for Module NAV_UKF

# Integration time
dt = 600 # [s]

# Process noise mapping
G = np.zeros((6, 3))
G[3:6, :] = np.eye(3)

# Acceleration noise level
sigma_a = 1e-6

# Continuous process noise covariance (spectral density)
Qc = sigma_a**2 * np.eye(3)

# Discrete process noise covariance
I3 = np.eye(3)
Q11 = (dt**3/3.0)*I3
Q12 = (dt**2/2.0)*I3
Q21 = Q12
Q22 = dt*I3

Q = sigma_a**2 * np.block([
    [Q11, Q12],
    [Q21, Q22]
])

NAV_UKF_par = {
    "dt"        : dt, # [s] Integration step. Update according to the SIM_dt
    "n_states"  : 6,  # Number of estimated states [pos[3], vel[3]]
    "G"         : G,  # Process noise mapping
    "Q"         : Q,  # Discrete process noise covariance
    "alpha"     : 1e-3,
    "beta"      : 2,
    "kappa"     : 0,

    # Initial values for the UKF
    "x_est_ini" : np.array([1.62228612e+08, 2.42235379e+07, 1.78903662e+07, -0.32883677, 27.86282899, 13.43536119]) + np.array([10, 10, 10, 0.01, 0.01, 0.01]), # Initial state
    "P_ini"     : np.eye(6), # Initial covariance
    "y_ini"     : 0, # Initial innovation
}
