# Level 2 Module DYN_EARTH
# Simulates the propagation of the Earth position relative to the SSB frame

import numpy as np
from Utils import spice
from .DYN_EARTH_par import DYN_EARTH_par

# Module output dictionary
def initialize(DYN_TIME_out):
    # Get initial position and orientation based on time
    time_UTC_ini = DYN_TIME_out["time_UTC"]
    EARTHpos_SSB_ini, EARTHvel_SSB_ini = spice.get_state("EARTH", time_UTC_ini)
    MOONpos_SSB_ini,  MOONvel_SSB_ini  = spice.get_state("MOON",  time_UTC_ini)
    # Update output
    DYN_EARTH_out = {
        "EARTHpos_SSB" : EARTHpos_SSB_ini,
        "EARTHvel_SSB" : EARTHvel_SSB_ini,
        "MOONpos_SSB"  : MOONpos_SSB_ini,
        "MOONvel_SSB"  : MOONvel_SSB_ini,
    }
    return DYN_EARTH_out

# Module main function
def outputs(t, DYN_out):
    EARTHpos_SSB, EARTHvel_SSB = spice.get_state("EARTH", DYN_out["DYN_TIME"]["time_UTC"])
    MOONpos_SSB, MOONvel_SSB   = spice.get_state("MOON",  DYN_out["DYN_TIME"]["time_UTC"])

    DYN_out["DYN_EARTH"]["EARTHpos_SSB"] = EARTHpos_SSB
    DYN_out["DYN_EARTH"]["EARTHvel_SSB"] = EARTHvel_SSB
    DYN_out["DYN_EARTH"]["MOONpos_SSB"]  = MOONpos_SSB
    DYN_out["DYN_EARTH"]["MOONvel_SSB"]  = MOONvel_SSB
 
    return DYN_out

# Module computation of derivatives to be integrated
def derivatives(t, DYN_out):
    return np.array([])

# Return integrated variables
def get_state(DYN_EARTH_out):
    return np.array([])

# Update integrated variables into the state dict
def set_state(DYN_EARTH_out, vec):
    return DYN_EARTH_out

