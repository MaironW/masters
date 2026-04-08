import numpy as np
from Utils import quaternions

# Temperature sensor orientation on BOF frame
CMBq_BOF = quaternions.rotvec2q(np.array([-np.pi/2, 0, 0]))

# Parameters for Module SEN_CMB
SEN_CMB_par = {
    # Temperature detector parameters
    "dt"         : 600,      # [s] Discretization of the sensor time - Also observation time
    "noise_mean" : 0,        # [K] Sensor noise mean
    "noise_std"  : 1e-6,     # [K] Sensor noise std
    "CMBq_BOF"   : CMBq_BOF, # Orientation of the STR relative to the Spacecraft body

    # Output initial values
    "SEN_CMBoutflg_ini" : 0,      # CMB (temperature) detector output flag
    "time_CMB_ini"      : 0,      # [s] CMB (temperature) detector time
    "T_dipole_ini"      : 2.7255, # [K] Measured temperatures
}
