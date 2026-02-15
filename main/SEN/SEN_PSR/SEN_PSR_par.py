import numpy as np

# Parameters for Module SEN_PSR

SEN_PSR_par = {
    "clock_bias" : 0, # [s] 

    # Output initial values
    "PSRoutflg_ini"    : 0,             # PSR (X-ray) detector output flag
    "time_PSR_ini"     : 0,             # PSR (X-ray) detector time
    "SCdt_SSB_mes_ini" : np.array([0]), # Measured time delay between the TOA of a pulse on the SC relative to the SSB (TDB)
}
