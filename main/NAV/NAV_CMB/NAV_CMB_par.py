import numpy as np
from Utils import quaternions

# Temperature sensors orientation on BOF frame (60 deg from spin axis, 120 deg apart, must be in sync with SEN_CMB)
CMB1q_BOF = quaternions.rotvec2q(np.pi/3 * np.array([0, 0, 1]))
CMB2q_BOF = quaternions.rotvec2q(np.pi/3 * np.array([0, -np.sin(np.deg2rad(120)), +np.cos(np.deg2rad(120))]))
CMB3q_BOF = quaternions.rotvec2q(np.pi/3 * np.array([0, -np.sin(np.deg2rad(240)), +np.cos(np.deg2rad(240))]))

# Parameters for Module NAV_CMB
NAV_CMB_par = {
    "sigma_T" : np.array([100e-6,
                          100e-6,
                          100e-6]), # [K] Temperature standard deviation per sensor
    "CMB1q_BOF"  : CMB1q_BOF, # Orientation of the CMB sensor 1 relative to the Spacecraft body
    "CMB2q_BOF"  : CMB2q_BOF, # Orientation of the CMB sensor 2 relative to the Spacecraft body
    "CMB3q_BOF"  : CMB3q_BOF, # Orientation of the CMB sensor 3 relative to the Spacecraft body
    "T_monopole" : 2.7255,    # [K] CMB average temperature

    # Output initial values
    "NAV_CMBoutflg_ini" : 0, # NAV_CMB Module output flag
    "z_ini"             : np.full(3, np.nan),
    "R_ini"             : np.full((3, 3), np.nan),
    "BOFq_SSB_ini"      : np.array([1,0,0,0]),
    "time_valid_ini"    : 0 # [s] initial time_CMB output for valid measurements
}
