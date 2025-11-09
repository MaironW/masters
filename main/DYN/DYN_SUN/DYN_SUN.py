# Level 2 Module DYN_SUN
# Simulates the propagation of the Sun position relative to the SSB frame

import numpy as np
from Utils import spice

from .DYN_SUN_par import DYN_SUN_par

# Module output dictionary
def initialize():
    DYN_SUN_out = {
        "SUNpos_SSB" : DYN_SUN_par["SUNpos_SSB_ini"],
        "SUNvel_SSB" : DYN_SUN_par["SUNvel_SSB_ini"],
        "SSBq_SUN"   : DYN_SUN_par["SSBq_SUN_ini"]
    }
    return DYN_SUN_out

# Module main function
def outputs(t, DYN_out):
    SUNpos_SSB, SUNvel_SSB = spice.get_state("SUN", DYN_out["DYN_TIME"]["time_UTC"])
    SSBq_SUN = spice.get_orientation("IAU_SUN", "J2000", DYN_out["DYN_TIME"]["time_UTC"])

    DYN_out["DYN_SUN"]["SUNpos_SSB"] = SUNpos_SSB
    DYN_out["DYN_SUN"]["SUNvel_SSB"] = SUNvel_SSB
    DYN_out["DYN_SUN"]["SSBq_SUN"]   = SSBq_SUN

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
