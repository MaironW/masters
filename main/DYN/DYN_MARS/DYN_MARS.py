# Level 2 Module DYN_MARS
# Simulates the propagation of the Mars position relative to the SSB frame

from Utils import spice

from .DYN_MARS_par import DYN_MARS_par

# Module output dictionary
def initialize():
    DYN_MARS_out = {
        "MARSpos_SSB"  : DYN_MARS_par["MARSpos_SSB_ini"],
        "MARSvel_SSB"  : DYN_MARS_par["MARSvel_SSB_ini"],
        "DEIMOSpos_SSB": DYN_MARS_par["DEIMOSpos_SSB_ini"],
        "DEIMOSvel_SSB": DYN_MARS_par["DEIMOSvel_SSB_ini"],
        "PHOBOSpos_SSB": DYN_MARS_par["PHOBOSpos_SSB_ini"],
        "PHOBOSvel_SSB": DYN_MARS_par["PHOBOSvel_SSB_ini"],
        "MCIq_MAR"     : DYN_MARS_par["MCIq_MAR_ini"]
    }
    return DYN_MARS_out

# Module main function
def outputs(t, DYN_out):
    MARSpos_SSB,   MARSvel_SSB   = spice.get_state("MARS",   DYN_out["DYN_TIME"]["time_UTC"])
    DEIMOSpos_SSB, DEIMOSvel_SSB = spice.get_state("DEIMOS", DYN_out["DYN_TIME"]["time_UTC"])
    PHOBOSpos_SSB, PHOBOSvel_SSB = spice.get_state("PHOBOS", DYN_out["DYN_TIME"]["time_UTC"])
    MCIq_MAR = spice.get_orientation("IAU_MARS", "J2000", DYN_out["DYN_TIME"]["time_UTC"])

    DYN_out["DYN_MARS"]["MARSpos_SSB"]   = MARSpos_SSB
    DYN_out["DYN_MARS"]["MARSvel_SSB"]   = MARSvel_SSB
    DYN_out["DYN_MARS"]["DEIMOSpos_SSB"] = DEIMOSpos_SSB
    DYN_out["DYN_MARS"]["PHOBOSpos_SSB"] = PHOBOSpos_SSB
    DYN_out["DYN_MARS"]["MCIq_MAR"]      = MCIq_MAR

    return DYN_out
