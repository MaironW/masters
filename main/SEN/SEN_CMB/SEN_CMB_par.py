import numpy as np
from Utils import quaternions

# Temperature sensors orientation on BOF frame (60 deg from spin axis, 120 deg apart)
# CSF = CMB Sensor Frame
z_angle = np.deg2rad(-30) # [rad]
y_angle = np.deg2rad(-60) # [rad]
qz = quaternions.rotvec2q([0,0,z_angle])
qy = quaternions.rotvec2q([0,y_angle,0])
CSF1q_BOF = quaternions.qprod(qy,qz)

z_angle = np.deg2rad(-150) # [rad]
y_angle = np.deg2rad(-60) # [rad]
qz = quaternions.rotvec2q([0,0,z_angle])
qy = quaternions.rotvec2q([0,y_angle,0])
CSF2q_BOF = quaternions.qprod(qy,qz)

z_angle = np.deg2rad(-270) # [rad]
y_angle = np.deg2rad(-60) # [rad]
qz = quaternions.rotvec2q([0,0,z_angle])
qy = quaternions.rotvec2q([0,y_angle,0])
CSF3q_BOF = quaternions.qprod(qy,qz)

# Parameters for Module SEN_CMB
SEN_CMB_par = {
    # Temperature detector parameters
    "dt"         : 600,       # [s] Discretization of the sensor time - Also observation time
    "noise_mean" : 0,         # [K] Sensor noise mean
    "noise_std"  : 100e-6,    # [K] Sensor noise std
    "CSF1q_BOF"  : CSF1q_BOF, # Orientation of the CMB sensor 1 relative to the Spacecraft body
    "CSF2q_BOF"  : CSF2q_BOF, # Orientation of the CMB sensor 2 relative to the Spacecraft body
    "CSF3q_BOF"  : CSF3q_BOF, # Orientation of the CMB sensor 3 relative to the Spacecraft body

    # Output initial values
    "SEN_CMBoutflg_ini" : 0,      # CMB (temperature) detector output flag
    "time_CMB_ini"      : 0,      # [s] CMB (temperature) detector time
    "T_dipole_ini"      : 2.7255, # [K] Measured temperatures
    "T_anisotropic_ini" : 0,      # [K] Measured anisotropies
}
