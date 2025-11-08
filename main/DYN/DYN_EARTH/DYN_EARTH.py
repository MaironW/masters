# Level 2 Module DYN_EARTH
# Simulates the propagation of the Earth position relative to the SSB frame

from Utils import spice

from .DYN_EARTH_par import DYN_EARTH_par

# Module output dictionary
def initialize():
    DYN_EARTH_out = {
        "EARTHpos_SSB" : DYN_EARTH_par["EARTHpos_SSB_ini"],
        "EARTHvel_SSB" : DYN_EARTH_par["EARTHvel_SSB_ini"],
        "MOONpos_SSB"  : DYN_EARTH_par["MOONpos_SSB_ini"],
        "ECIq_TER"     : DYN_EARTH_par["ECIq_TER_ini"],
        "TERq_ECI"     : DYN_EARTH_par["TERq_ECI_ini"],
    }
    return DYN_EARTH_out

# Module main function
def outputs(t, DYN_out):
    EARTHpos_SSB, EARTHvel_SSB = spice.get_state("EARTH", DYN_out["DYN_TIME"]["time_UTC"])
    MOONpos_SSB, MOONvel_SSB   = spice.get_state("MOON", DYN_out["DYN_TIME"]["time_UTC"])
    ECIq_TER = spice.get_orientation("ITRF93", "J2000", DYN_out["DYN_TIME"]["time_UTC"])
    TERq_ECI = spice.get_orientation("J2000", "ITRF93", DYN_out["DYN_TIME"]["time_UTC"])

    DYN_out["DYN_EARTH"]["EARTHpos_SSB"] = EARTHpos_SSB
    DYN_out["DYN_EARTH"]["EARTHvel_SSB"] = EARTHvel_SSB
    DYN_out["DYN_EARTH"]["MOONpos_SSB"]  = MOONpos_SSB
    DYN_out["DYN_EARTH"]["ECIq_TER"]     = ECIq_TER
    DYN_out["DYN_EARTH"]["TERq_ECI"]     = TERq_ECI
 
    return DYN_out
