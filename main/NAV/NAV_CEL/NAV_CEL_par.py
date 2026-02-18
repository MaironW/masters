import numpy as np
from Utils.constants import CONSTANTS_par

# Pre-computations
los_angle_min = 10*CONSTANTS_par["deg2rad_cst"]  # [rad] Minimum angle between body and star to be considered valid
los_angle_max = 120*CONSTANTS_par["deg2rad_cst"] # [rad] Maximum angle between body and star to be considered valid
cos_los_angle_min = np.cos(los_angle_max)
cos_los_angle_max = np.cos(los_angle_min)

# Parameters for Module NAV_CEL
NAV_CEL_par = {
    # Celestial Navigation Parameters
    "cos_los_angle_min" : cos_los_angle_min, # Minimum cosine for angle between body and star to be considered valid
    "cos_los_angle_max" : cos_los_angle_max, # Maximum cosine for angle between body and star to be considered valid
    "n_bodies"          : 6,                 # Number of celestial bodies to be evaluated as reference
    # Outputs initial values
    "NAV_CELoutflg_ini"       : 0,                      # NAV_CEL Module output flag
    "z_ini"                   : np.full(2, np.nan),     # Measurement vector
    "R_ini"                   : np.full((2,2), np.nan), # Measurement covariance
    "sigma_angle"             : 1e-3,                   # [rad] STR noise std in [X,Y,Z] directions
    "BODYangles_mes_ini"      : np.full(1, np.nan),     # [rad] LOS angle between the body and one star
    "BODYsel_STARdir_mes_ini" : np.full(3, np.nan),     # LOS direction between the SC and one star
}
