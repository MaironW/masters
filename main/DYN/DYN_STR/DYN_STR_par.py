import numpy as np

# Parameters for Module DYN_STR
DYN_STR_par = {
    "STARSdir_SSB_ini" : np.array([0,0,0]),  # Initial direction of the stars in the SSB frame at the initial epoch
    "STARSid_ini"      : np.array([np.nan]), # Initial direction of the stars in the SSB frame at the initial epoch
    "magnitude_min"    : 6,                  # DYN_STR will only select stars with higher brightness than this value
    "num_stars_max"    : 250,                # Maximum number of stars to load 
}
