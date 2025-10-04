# Level 2 Module DYN_MARS
# Simulates the propagation of the Mars position relative to the SSB frame

from Utils import spice

from .DYN_MARS_par import DYN_MARS_par

# Module output dictionary
DYN_MARS_out = {
    "MARSpos_SSB"  : DYN_MARS_par["MARSpos_SSB_ini"],
    "MARSvel_SSB"  : DYN_MARS_par["MARSvel_SSB_ini"],
    "DEIMOSpos_SSB": DYN_MARS_par["DEIMOSpos_SSB_ini"],
    "DEIMOSvel_SSB": DYN_MARS_par["DEIMOSvel_SSB_ini"],
    "PHOBOSpos_SSB": DYN_MARS_par["PHOBOSpos_SSB_ini"],
    "PHOBOSvel_SSB": DYN_MARS_par["PHOBOSvel_SSB_ini"],
}

# Module main function
def run(DYN_TIME_out):
    MARSpos_SSB,   MARSvel_SSB   = spice.get_state("MARS",   DYN_TIME_out["time_UTC"])
    DEIMOSpos_SSB, DEIMOSvel_SSB = spice.get_state("DEIMOS", DYN_TIME_out["time_UTC"])
    PHOBOSpos_SSB, PHOBOSvel_SSB = spice.get_state("PHOBOS", DYN_TIME_out["time_UTC"])

    DYN_MARS_out["MARSpos_SSB"]   = MARSpos_SSB
    DYN_MARS_out["MARSvel_SSB"]   = MARSvel_SSB
    DYN_MARS_out["DEIMOSpos_SSB"] = DEIMOSpos_SSB
    DYN_MARS_out["PHOBOSpos_SSB"] = PHOBOSpos_SSB

    return dict(DYN_MARS_out)
