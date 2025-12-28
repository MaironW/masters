import numpy as np
from Utils.constants import CONSTANTS_par

# Parameters for Module DYN_TRA

DYN_TRA_par = {
    "ref_elements"  : "kep", # Use 'kep' elements or 'rvi' state vectors to initialize DYN_TRA
    "BODY_ini"      : "SUN", # Select with respect to which body the orbital elements are defined
    "sma_ini"       : 204687945.33840084,                              # [km] Semi-major axis
    "ecc_ini"       : 0.2578141373442776,                              # Eccentricity
    "incl_ini"      : 25.896289319887888*CONSTANTS_par["deg2rad_cst"], # [rad] Inclination
    "raan_ini"      : 355.4062892154075*CONSTANTS_par["deg2rad_cst"],  # [rad] Right ascension of the ascending node
    "argp_ini"      : 321.24872768327856*CONSTANTS_par["deg2rad_cst"], # [rad] Argument of perigee
    "tano_ini"      : 5.647168512920133*CONSTANTS_par["deg2rad_cst"],  # [rad] True anomaly
    "SCpos_SSB_ini" : np.array([0,0,0]),                               # [km]   Initial position relative to SSB frame
    "SCvel_SSB_ini" : np.array([0,0,0]),                               # [km/s] Initial velocity relative to SSB frame
}
