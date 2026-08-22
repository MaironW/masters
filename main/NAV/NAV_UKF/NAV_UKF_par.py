import numpy as np

# Parameters for Module NAV_UKF

# Integration time
dt = 60 # [s] 1 min

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

# Initial conditions based on the same values defined on DYN_TRA_par
SCpos_SSB_true_ini = np.array([-3.10735179e+07, 1.41664640e+08, 5.83463713e+07]) # [km]   Initial position relative to SSB frame
SCvel_SSB_true_ini = np.array([-31.44040794, -1.35567615, -1.45138134])          # [km/s] Initial velocity relative to SSB frame

SCpos_SSB_err_ini = np.array([10,10,10]) # [km]
SCvel_SSB_err_ini = np.array([0.01,0.01,0.01]) # [km/s]

SCpos_SSB_est_ini = SCpos_SSB_true_ini + SCpos_SSB_err_ini
SCvel_SSB_est_ini = SCvel_SSB_true_ini + SCvel_SSB_err_ini

x_est_ini = np.concatenate((SCpos_SSB_est_ini, SCvel_SSB_est_ini))

NAV_UKF_par = {
    "dt"        : dt, # [s] Integration step. Update according to the SIM_dt
    "n_states"  : 6,  # Number of estimated states [pos[3], vel[3]]
    "G"         : G,  # Process noise mapping
    "Q"         : Q,  # Discrete process noise covariance
    "alpha"     : 1e-3,
    "beta"      : 2,
    "kappa"     : 0,

    # Initial values for the UKF
    "x_est_ini" : x_est_ini,
    "P_ini"     : np.eye(6), # Initial covariance
    "y_ini"     : 0, # Initial innovation
}
