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

    # X-ray properties
    "background_photon_flux" : 50,    # [photon/m^2/s] Background photon flux
    "pulsar_photon_flux"     : 1.93,  # [photon/m^2/s] Pulsar photon flux
    "pulse_width"            : 0.3e-3, # [s]

    # Output initial values
    "PULSARSdir_SSB_ini" : np.array([0,0,0]), # Initial direction of the pulsars in the SSB frame at the initial epoch
    "phase_SSB_ini"      : np.array([0]), # Initial phase for the pulsars measured in the SSB
    "phase_SC_ini"       : np.array([0]), # Initial phase for the pulsars measured in the SC
    "SCdt_SSB_ini"       : np.array([0]), # Time delay between the TOA of a pulse on the SC relative to the SSB (TDB)
}
