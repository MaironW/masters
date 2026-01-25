
import numpy as np
from Utils.constants import CONSTANTS_par

# Parameters for Module NAV_CEL
NAV_CEL_par = {
    # Celestial Navigation Parameters
    "los_angle_min" : 1*CONSTANTS_par["deg2rad_cst"], # [rad] Minimum angle between body and star to be considered valid
    # Outputs initial values
    "NAV_CELoutflg_ini"  : 0,                 # NAV_CEL Module output flag
    "BODYangles_mes_ini" : np.full(3, np.nan) # [rad] LOS angles between the body and three stars
}
