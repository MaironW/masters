from Utils import misc
import numpy as np

# Define general constants which might be used through the entire code

CONSTANTS_par = {
    # Convertion factors
    "sec2days_cst" : 1/(3600*24),       # Seconds to days
    "deg2rad_cst"  :  0.017453293,      # Degrees to radians
    "rad2deg_cst"  : 57.295779513,      # Radians to degrees
    "day2sec_cst"  : 86400,             # Days to seconds
    "pc2m_cst"     : 3.08567758128e+16, # Parsec to meters
    "m2cm_cst"     : 100,               # Meters to centimeters
    # Mass properties
    "SCmass_cst"    : 1000.0,    # [kg] Spacecraft mass
    "SUNmass_cst"   : 1.9885e30, # [kg] Sun mass
    "EARTHmass_cst" : 5.9722e24, # [kg] Earth mass
    "MARSmass_cst"  : 6.4171e23, # [kg] Mars mass
    # Body sizes
    "SUNradius_cst"   : 696340.0, # [km] Sun radius
    "EARTHradius_cst" :   6371.0, # [km] Earth radius
    "MARSradius_cst"  :   3389.5, # [km] Mars radius
    # Universal constants
    "gravitational_cst" : 6.67430e-20, # [km^3/(kg s^2)]
    "light_speed_cst"   : 299792.458,  # [km/s] Light speed
    # Astronomical Unit
    "AU_cst" : 188766877, # [km]
    # Time
    "MJD2000epoch_relMJD_TDB_days_cst" : 51544.5, # MJD days to J2000 days
    # CMB temperature
    "T_CMBR_cst" : 2.7255, # [K] CMBR average temperature
}

# Compute derived constants
CONSTANTS_par["mu_SUN_cst"]   = CONSTANTS_par["gravitational_cst"]*CONSTANTS_par["SUNmass_cst"]   # [km^3/s^2] Sun gravitational parameter
CONSTANTS_par["mu_EARTH_cst"] = CONSTANTS_par["gravitational_cst"]*CONSTANTS_par["EARTHmass_cst"] # [km^3/s^2] Earth gravitational parameter
CONSTANTS_par["mu_MARS_cst"]  = CONSTANTS_par["gravitational_cst"]*CONSTANTS_par["MARSmass_cst"]  # [km^3/s^2] Mars gravitational parameter

# Velocity of the Solar System with respect to the CMBR on the SSB frame
SSBvel_GAL_norm  = 371 # [km/s] Absolute velocity of the Solar System
SSBvel_GAL_dir_l = 263.85*CONSTANTS_par["deg2rad_cst"] # [rad] Galactic longitude of the SSB velocity
SSBvel_GAL_dir_b = 48.25*CONSTANTS_par["deg2rad_cst"]  # [rad] Galactic latitude of the SSB velocity
SSBvel_GAL_dir   = misc.latlon2dir(SSBvel_GAL_dir_b, SSBvel_GAL_dir_l)
SSBvel_CMB       = misc.GALtoSSB(SSBvel_GAL_dir)
CONSTANTS_par["SSBvel_CMB_cst"] = SSBvel_CMB
