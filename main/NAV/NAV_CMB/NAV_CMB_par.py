import numpy as np
from Utils import quaternions

# Temperature sensors orientation on BOF frame (60 deg from spin axis, 120 deg apart, must be in sync with SEN_CMB)
z_angle = np.deg2rad(-150) # [rad]
y_angle = np.deg2rad(-60) # [rad]
qz = quaternions.rotvec2q([0,0,z_angle])
qy = quaternions.rotvec2q([0,y_angle,0])
CSF1q_BOF = quaternions.qprod(qy,qz)

z_angle = np.deg2rad(-30) # [rad]
y_angle = np.deg2rad(-60) # [rad]
qz = quaternions.rotvec2q([0,0,z_angle])
qy = quaternions.rotvec2q([0,y_angle,0])
CSF2q_BOF = quaternions.qprod(qy,qz)

z_angle = np.deg2rad(-270) # [rad]
y_angle = np.deg2rad(-60) # [rad]
qz = quaternions.rotvec2q([0,0,z_angle])
qy = quaternions.rotvec2q([0,y_angle,0])
CSF3q_BOF = quaternions.qprod(qy,qz)

# Parameters for Module NAV_CMB
NAV_CMB_par = {
    "sigma_T" : np.array([200e-6,
                          200e-6,
                          200e-6]), # [K] Temperature standard deviation per sensor
    "CSF1q_BOF"  : CSF1q_BOF, # Orientation of the CMB sensor 1 relative to the Spacecraft body
    "CSF2q_BOF"  : CSF2q_BOF, # Orientation of the CMB sensor 2 relative to the Spacecraft body
    "CSF3q_BOF"  : CSF3q_BOF, # Orientation of the CMB sensor 3 relative to the Spacecraft body
    "T_monopole" : 2.7255,    # [K] CMB average temperature

    # Output initial values
    "NAV_CMBoutflg_ini" : 0, # NAV_CMB Module output flag
    "z_ini"             : np.full(3, np.nan),
    "R_ini"             : np.full((3, 3), np.nan),
    "BOFq_SSB_ini"      : np.array([1,0,0,0]),
    "time_valid_ini"    : 0 # [s] initial time_CMB output for valid measurements
}
