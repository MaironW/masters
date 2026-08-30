import numpy as np

# Parameters for Module NAV_PSR
NAV_PSR_par = {
    # [s] TOA standard deviation per pulsar, tunned to reduce overconfidence
    # Base values computed in SEN_PSR
    "sigma_TOA" : np.array([1.39812310e-06,
                            1.71242967e-06,
                            1.75909018e-06,
                            9.26570994e-06,
                            1.01698121e-05,
                            1.60063706e-05,
                            7.08752914e-05,
                            8.68895977e-05,
                            4.58684526e-05,
                            8.42690951e-05]),

    # Limit the number of pulsars used in the simulation
    # It aways select the first pulsars (should be between 3 and 10 for valid navigation)
    "num_pulsars_max" : 3,

    # Output initial values
    "NAV_PSRoutflg_ini"  : 0, # NAV_PSR Module output flag
    "SSBpos_SUN_ref_ini" : np.array([0,0,0]),
    "time_valid_ini"     : -np.inf # [s] initial time_PSR output for valid measurements
}
