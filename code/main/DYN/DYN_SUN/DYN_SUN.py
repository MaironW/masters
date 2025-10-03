# Level 2 Module DYN_SUN
# Simulates the propagation of the Sun position relative to the SSB frame

from .DYN_SUN_par import DYN_SUN_par

# Module output dictionary
DYN_SUN_out = {
    "SUNpos_SSB" : DYN_SUN_par["SUNpos_SSB_ini"]
}

# Module main function
def run():
    DYN_SUN_out["SUNpos_SSB"] = DYN_SUN_par["SUNpos_SSB_ini"]
    return dict(DYN_SUN_out)
