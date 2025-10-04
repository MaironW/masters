# Level 2 Module DYN_SUN
# Simulates the propagation of the Sun position relative to the SSB frame

from Utils import spice

from .DYN_SUN_par import DYN_SUN_par

# Module output dictionary
DYN_SUN_out = {
    "SUNpos_SSB" : DYN_SUN_par["SUNpos_SSB_ini"],
    "SUNvel_SSB" : DYN_SUN_par["SUNvel_SSB_ini"]
}

# Module main function
def run(DYN_TIME_out):
    SUNpos_SSB, SUNvel_SSB = spice.get_state("SUN", DYN_TIME_out["time_UTC"])

    DYN_SUN_out["SUNpos_SSB"] = SUNpos_SSB
    DYN_SUN_out["SUNvel_SSB"] = SUNvel_SSB

    return dict(DYN_SUN_out)
