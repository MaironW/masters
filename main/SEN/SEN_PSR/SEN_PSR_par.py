import numpy as np

# Parameters for Module SEN_PSR
SEN_PSR_par = {
    # X-ray detector parameters
    "dt"            : 600,  # [s]   Discretization of the sensor time - Also observation time
    "detector_area" : 0.19, # [m^2] Detector area
    "t_bias"        : 0,    # [s]   Detector (spacecraft) clock bias
    "SNR_max"       : 1000, # Maximuim SNR value

    # Output initial values
    "SEN_PSRoutflg_ini"     : 0,                # PSR (X-ray) detector output flag
    "time_PSR_ini"          : 0,                # PSR (X-ray) detector time
    "OBTdt_TDB_mes_ini"     : np.array([0]),    # [s] Measured time delay between the TOA of a pulse on the SC (OBT) relative to the SSB (TDB)
    "covariance_ini"        : np.array([0]),    # [s] Sensor time delay covariance
    "PULSARSdir_SC_mes_ini" : np.array([0]),    # Measured pulsar direction from the SC
    "PULSARSid_mes_ini"     : np.array([None]), # Dummy, to be filled at initialization
}
