# Define general constants which might be used through the entire code

CONSTANTS_par = {
    # Convertion factors
    "sec2days_cst" : 1/(3600*24),  # Seconds to days
    "deg2rad_cst"  :  0.017453293, # Degrees to radians
    "rad2deg_cst"  : 57.295779513, # Radians to degrees
    "day2sec_cst"  : 86400,        # Days to seconds
    # Mass properties
    "SCmass_cst"    : 1000.0,    # [kg] Spacecraft mass
    "SUNmass_cst"   : 1.9885e30, # [kg] Sun mmass
    "EARTHmass_cst" : 5.9722e24, # [kg] Earth mass
    "MARSmass_cst"  : 6.4171e23, # [kg] Mars mass
    # Body sizes
    "SUNradius_cst"   : 696340.0, # [km] Sun radius
    "EARTHradius_cst" :   6371.0, # [km] Earth radius
    "MARSradius_cst"  :   3389.5, # [km] Mars radius
    # Universal constants
    "gravitational_cst" : 6.67430e-20, # [km^3/(kg s^2)]
    "light_speed_cst" : 299792.458,    # [km/s] Light speed
    # Astronomical Unit
    "AU_cst" : 188766877, # [km]
    # Time
    "MJD2000epoch_relMJD_TDB_days" : 51544.5, # MJD days to J2000 days
}

# Compute derived constants
CONSTANTS_par["mu_SUN_cst"]   = CONSTANTS_par["gravitational_cst"]*CONSTANTS_par["SUNmass_cst"]   # [km^3/s^2] Sun gravitational parameter
CONSTANTS_par["mu_EARTH_cst"] = CONSTANTS_par["gravitational_cst"]*CONSTANTS_par["EARTHmass_cst"] # [km^3/s^2] Earth gravitational parameter
CONSTANTS_par["mu_MARS_cst"]  = CONSTANTS_par["gravitational_cst"]*CONSTANTS_par["MARSmass_cst"]  # [km^3/s^2] Mars gravitational parameter