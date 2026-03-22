import numpy as np

# Parameters for Module NAV_PSR
NAV_PSR_par = {
    "sigma_TOA" : np.array([1.04913991e-06,
                            5.68710956e-06,
                            4.35825866e-05,
                            1.06742387e-06,
                            1.27842941e-05,
                            3.94597966e-04,
                            2.26918341e-04,
                            5.93848437e-05,
                            5.34729613e-05,
                            3.37885772e-07]), # [s] TOA standard deviation per pulsar

    # Output initial values
    "NAV_PSRoutflg_ini"  : 0, # NAV_PSR Module output flag
    "SSBpos_SUN_ref_ini" : np.array([0,0,0])
}
