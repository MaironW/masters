import numpy as np
from Utils import quaternions

# Temperature sensors orientation on BOF frame
CMB1q_BOF = quaternions.rotvec2q(np.pi/6 * np.array([0, 0, 1]))
CMB2q_BOF = quaternions.rotvec2q(np.pi/6 * np.array([0, -np.sin(np.deg2rad(120)), +np.cos(np.deg2rad(120))]))
CMB3q_BOF = quaternions.rotvec2q(np.pi/6 * np.array([0, -np.sin(np.deg2rad(240)), +np.cos(np.deg2rad(240))]))

# Parameters for Module SEN_CMB
SEN_CMB_par = {
    # Temperature detector parameters
    "dt"         : 600,       # [s] Discretization of the sensor time - Also observation time
    "noise_mean" : 0,         # [K] Sensor noise mean
    "noise_std"  : 1e-5,      # [K] Sensor noise std
    "CMB1q_BOF"  : CMB1q_BOF, # Orientation of the CMB sensor 1 relative to the Spacecraft body
    "CMB2q_BOF"  : CMB2q_BOF, # Orientation of the CMB sensor 2 relative to the Spacecraft body
    "CMB3q_BOF"  : CMB3q_BOF, # Orientation of the CMB sensor 3 relative to the Spacecraft body

    # Output initial values
    "SEN_CMBoutflg_ini" : 0,      # CMB (temperature) detector output flag
    "time_CMB_ini"      : 0,      # [s] CMB (temperature) detector time
    "T_dipole_ini"      : 2.7255, # [K] Measured temperatures
}
