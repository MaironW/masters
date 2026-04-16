import numpy as np

# Parameters for Module NAV_EKF

NAV_EKF_par = {
    "n_states"  : 6, # Number of estimated states [pos[3], vel[3]]

    # Initial values for the EKF
    "x_est_ini" : np.array([0,0,0,0,0,0]),
    "P_ini"     : np.ones((6,6)),
    "y_inn_ini" : np.array([0,0,0,0,0,0]),
}
