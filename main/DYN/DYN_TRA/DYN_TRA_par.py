import numpy as np
from Utils.constants import CONSTANTS_par

# Parameters for Module DYN_TRA

DYN_TRA_par = {
    "sma_ini"  : CONSTANTS_par["EARTHradius_cst"] + 1000, # [km] Semi-major axis
    "ecc_ini"  : 0,                                       # Eccentricity
    "incl_ini" : CONSTANTS_par["deg2rad_cst"]*99.5,       # [rad] Inclination
    "raan_ini" : 0,                                       # [rad] Right ascension of the ascending node
    "argp_ini" : 0,                                       # [rad] Argument of perigee
    "tano_ini" : 0,                                       # [rad] True anomaly
}