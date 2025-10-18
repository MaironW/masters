# Define general constants which might be used through the entire code

CONSTANTS_par = {
    # Convertion factors
    "sec2days_cst"  : 1/(3600*24),
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
    "gravitational_cst" : 6.67430e-20 # [km^3/(kg s^2)]
}