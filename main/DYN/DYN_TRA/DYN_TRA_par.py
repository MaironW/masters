import numpy as np
from Utils.constants import CONSTANTS_par

# Parameters for Module DYN_TRA

DYN_TRA_par = {
    "ref_elements"  : "rvi", # Use 'kep' elements or 'rvi' state vectors to initialize DYN_TRA
    "BODY_ini"      : "SUN", # Select with respect to which body the orbital elements are defined
    "sma_ini"       : None, # [km] Semi-major axis
    "ecc_ini"       : None, # Eccentricity
    "incl_ini"      : None, # [rad] Inclination
    "raan_ini"      : None, # [rad] Right ascension of the ascending node
    "argp_ini"      : None, # [rad] Argument of perigee
    "tano_ini"      : None, # [rad] True anomaly
    # MSL orbit initial conditions, at "2012-01-01 T00:00:00"
    "SCpos_SSB_ini" : np.array([-3.10735179e+07, 1.41664640e+08, 5.83463713e+07]), # [km]   Initial position relative to SSB frame
    "SCvel_SSB_ini" : np.array([-31.44040794, -1.35567615, -1.45138134]),          # [km/s] Initial velocity relative to SSB frame
    "SCacc_SSB_ini" : np.array([0, 0, 0]), # [km/s^2] Initial acceleration relative to SSB frame
}
