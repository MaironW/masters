# Level 2 Module DYN_EARTH
# Simulates the propagation of the Earth position relative to the SSB frame

from Utils import spice

from .DYN_EARTH_par import DYN_EARTH_par

# Module output dictionary
DYN_EARTH_out = {
    "EARTHpos_SSB" : DYN_EARTH_par["EARTHpos_SSB_ini"],
    "EARTHvel_SSB" : DYN_EARTH_par["EARTHvel_SSB_ini"],
    "MOONpos_SSB"  : DYN_EARTH_par["MOONpos_SSB_ini"],
    "ECIq_TER"     : DYN_EARTH_par["ECIq_TER_ini"],
    "TERq_ECI"     : DYN_EARTH_par["TERq_ECI_ini"],
}

# Module main function
def run(DYN_TIME_out):
    EARTHpos_SSB, EARTHvel_SSB = spice.get_state("EARTH", DYN_TIME_out["time_UTC"])
    MOONpos_SSB, MOONvel_SSB   = spice.get_state("MOON", DYN_TIME_out["time_UTC"])
    ECIq_TER = spice.get_orientation("IAU_EARTH", "J2000", DYN_TIME_out["time_UTC"])
    TERq_ECI = spice.get_orientation("J2000", "IAU_EARTH", DYN_TIME_out["time_UTC"])

    DYN_EARTH_out["EARTHpos_SSB"] = EARTHpos_SSB
    DYN_EARTH_out["EARTHvel_SSB"] = EARTHvel_SSB
    DYN_EARTH_out["MOONpos_SSB"]  = MOONpos_SSB
    DYN_EARTH_out["ECIq_TER"]     = ECIq_TER
    DYN_EARTH_out["TERq_ECI"]     = TERq_ECI
 
    return dict(DYN_EARTH_out)
