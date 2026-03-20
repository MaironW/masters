import numpy as np

# Parameters for Module DYN_PSR
DYN_PSR_par = {
    # Select pulsars with the criteria:
    # - millisecond pulsar (high frequency)
    # - high enery pulsars
    # - non-binary pulsars
    "selection_criteria" : "f0 > 50  && \
                            TYPE(HE) && \
                            !TYPE(BINARY) &&\
                            exist(PEPOCH) && \
                            exist(F0) && \
                            exist(F1) && \
                            exist(RAJD) && \
                            exist(DECJD)",

    # Output initial values
    "PULSARSdir_SSB_ini" : np.array([0,0,0]), # Initial direction of the pulsars in the SSB frame at the initial epoch
    "roemer_delay_ini"   : np.array([0]), # [s] Time delay due to Doppler and annual parallax effect
    "shapiro_delay_ini"  : np.array([0]), # [s] Time delay due to gravity distortion of light
    "OBTdt_TDB_ini"      : np.array([0]), # [s] Time delay between the TOA of a pulse on the SC (OBT) relative to the SSB (TDB)
}
