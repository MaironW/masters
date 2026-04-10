import numpy as np
from Utils import quaternions

# Temperature sensors orientation on BOF frame
CMB1q_BOF = quaternions.rotvec2q(np.pi/6 * np.array([0, 0, 1]))
CMB2q_BOF = quaternions.rotvec2q(np.pi/6 * np.array([0, -np.sin(np.deg2rad(120)), +np.cos(np.deg2rad(120))]))
CMB3q_BOF = quaternions.rotvec2q(np.pi/6 * np.array([0, -np.sin(np.deg2rad(240)), +np.cos(np.deg2rad(240))]))

# Parameters for Module NAV_CMB
NAV_CMB_par = {
    "sigma_T" : np.array([1.04913991e-06,
                            5.68710956e-06,
                            4.35825866e-05]), # [s] Temperature standard deviation per sensor
    "CMB1q_BOF"  : CMB1q_BOF, # Orientation of the CMB sensor 1 relative to the Spacecraft body
    "CMB2q_BOF"  : CMB2q_BOF, # Orientation of the CMB sensor 2 relative to the Spacecraft body
    "CMB3q_BOF"  : CMB3q_BOF, # Orientation of the CMB sensor 3 relative to the Spacecraft body

    # Output initial values
    "NAV_CMBoutflg_ini"  : 0, # NAV_CMB Module output flag
    "z_ini"              : np.full(3, np.nan),
    "R_ini"              : np.full((3, 3), np.nan),
}
