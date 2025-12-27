# Level 2 Module DYN_SUN
# Simulates the propagation of the Sun position relative to the SSB frame

import numpy as np
from Utils import spice

from .DYN_SUN_par import DYN_SUN_par

# Module output dictionary
def initialize(DYN_TIME_out):
    # Get initial position and orientation based on time
    time_UTC_ini = DYN_TIME_out["time_UTC"]
    SUNpos_SSB_ini, SUNvel_SSB_ini = spice.get_state("SUN", time_UTC_ini)
    # Update output
    DYN_SUN_out = {
        "SUNpos_SSB" : SUNpos_SSB_ini,
        "SUNvel_SSB" : SUNvel_SSB_ini,
    }
    return DYN_SUN_out

# Module main function
def outputs(t, DYN_out):
    SUNpos_SSB, SUNvel_SSB = spice.get_state("SUN", DYN_out["DYN_TIME"]["time_UTC"])

    DYN_out["DYN_SUN"]["SUNpos_SSB"] = SUNpos_SSB
    DYN_out["DYN_SUN"]["SUNvel_SSB"] = SUNvel_SSB

    return DYN_out

# Module computation of derivatives to be integrated
def derivatives(t, DYN_out):
    return np.array([])

# Return integrated variables
def get_state(DYN_SUN_out):
    return np.array([])

# Update integrated variables into the state dict
def set_state(DYN_SUN_out, vec):
    return DYN_SUN_out
