import numpy as np

# Parameters for Module NAV_PSR
NAV_PSR_par = {
    "sigma_TOA" : np.array([0.01058292, 0.06837116, 0.01218063, 0.00420778]), # [s] TOA standard deviation

    # Output initial values
    "NAV_PSRoutflg_ini"  : 0, # NAV_PSR Module output flag
    "SSBpos_SUN_ref_ini" : np.array([0,0,0])
}
