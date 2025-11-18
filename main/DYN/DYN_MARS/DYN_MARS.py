# Level 2 Module DYN_MARS
# Simulates the propagation of the Mars position relative to the SSB frame

import numpy as np
from Utils import spice

from .DYN_MARS_par import DYN_MARS_par

# Module output dictionary
def initialize(DYN_TIME_out):
    # Get initial position and orientation based on time
    time_UTC_ini = DYN_TIME_out["time_UTC"]
    MARSpos_SSB_ini,   MARSvel_SSB_ini   = spice.get_state("MARS",   time_UTC_ini)
    DEIMOSpos_SSB_ini, DEIMOSvel_SSB_ini = spice.get_state("DEIMOS", time_UTC_ini)
    PHOBOSpos_SSB_ini, PHOBOSvel_SSB_ini = spice.get_state("PHOBOS", time_UTC_ini)
    MCIq_MAR_ini = spice.get_orientation("IAU_MARS", "J2000", time_UTC_ini)
    MARq_MCI_ini = spice.get_orientation("J2000", "IAU_MARS", time_UTC_ini)
    # Update output
    DYN_MARS_out = {
        "MARSpos_SSB"  : MARSpos_SSB_ini,
        "MARSvel_SSB"  : MARSvel_SSB_ini,
        "DEIMOSpos_SSB": DEIMOSpos_SSB_ini,
        "DEIMOSvel_SSB": DEIMOSvel_SSB_ini,
        "PHOBOSpos_SSB": PHOBOSpos_SSB_ini,
        "PHOBOSvel_SSB": PHOBOSvel_SSB_ini,
        "MCIq_MAR"     : MCIq_MAR_ini,
        "MARq_MCI"     : MARq_MCI_ini
    }
    return DYN_MARS_out

# Module main function
def outputs(t, DYN_out):
    MARSpos_SSB,   MARSvel_SSB   = spice.get_state("MARS",   DYN_out["DYN_TIME"]["time_UTC"])
    DEIMOSpos_SSB, DEIMOSvel_SSB = spice.get_state("DEIMOS", DYN_out["DYN_TIME"]["time_UTC"])
    PHOBOSpos_SSB, PHOBOSvel_SSB = spice.get_state("PHOBOS", DYN_out["DYN_TIME"]["time_UTC"])
    MCIq_MAR = spice.get_orientation("IAU_MARS", "J2000", DYN_out["DYN_TIME"]["time_UTC"])
    MARq_MCI = spice.get_orientation("J2000", "IAU_MARS", DYN_out["DYN_TIME"]["time_UTC"])

    DYN_out["DYN_MARS"]["MARSpos_SSB"]   = MARSpos_SSB
    DYN_out["DYN_MARS"]["MARSvel_SSB"]   = MARSvel_SSB
    DYN_out["DYN_MARS"]["DEIMOSpos_SSB"] = DEIMOSpos_SSB
    DYN_out["DYN_MARS"]["DEIMOSvel_SSB"] = DEIMOSvel_SSB
    DYN_out["DYN_MARS"]["PHOBOSpos_SSB"] = PHOBOSpos_SSB
    DYN_out["DYN_MARS"]["PHOBOSvel_SSB"] = PHOBOSvel_SSB
    DYN_out["DYN_MARS"]["MCIq_MAR"]      = MCIq_MAR
    DYN_out["DYN_MARS"]["MARq_MCI"]      = MARq_MCI

    return DYN_out

# Module computation of derivatives to be integrated
def derivatives(t, DYN_out):
    return np.array([])

# Return integrated variables
def get_state(DYN_MARS_out):
    return np.array([])

# Update integrated variables into the state dict
def set_state(DYN_MARS_out, vec):
    return DYN_MARS_out
