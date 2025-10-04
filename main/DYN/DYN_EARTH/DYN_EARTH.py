# Level 2 Module DYN_EARTH
# Simulates the propagation of the Earth position relative to the SSB frame

from Utils import spice

from .DYN_EARTH_par import DYN_EARTH_par

# Module output dictionary
DYN_EARTH_out = {
    "EARTHpos_SSB" : DYN_EARTH_par["EARTHpos_SSB_ini"],
    "EARTHvel_SSB" : DYN_EARTH_par["EARTHvel_SSB_ini"],
    "MOONpos_SSB"  : DYN_EARTH_par["MOONpos_SSB_ini"],
}

# Module main function
def run(DYN_TIME_out):
    EARTHpos_SSB, EARTHvel_SSB = spice.get_state("EARTH", DYN_TIME_out["time_UTC"])
    MOONpos_SSB, MOONvel_SSB   = spice.get_state("MOON", DYN_TIME_out["time_UTC"])

    DYN_EARTH_out["EARTHpos_SSB"] = EARTHpos_SSB
    DYN_EARTH_out["EARTHvel_SSB"] = EARTHvel_SSB
    DYN_EARTH_out["MOONpos_SSB"]  = MOONpos_SSB

    return dict(DYN_EARTH_out)
