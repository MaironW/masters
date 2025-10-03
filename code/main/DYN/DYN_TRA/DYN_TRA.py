# Level 2 Module DYN_TRA
# Simulates the propagation of the spacecraft orbital states for the simulation

from .DYN_TRA_par import DYN_TRA_par

# Module output dictionary
DYN_TRA_out = {
    "SCpos_SSB" : DYN_TRA_par["SCpos_SSB_ini"],
    "SCvel_SSB" : DYN_TRA_par["SCvel_SSB_ini"]
}

# Module main function
def run():
    DYN_TRA_out["SCpos_SSB"] = DYN_TRA_par["SCpos_SSB_ini"]
    DYN_TRA_out["SCvel_SSB"] = DYN_TRA_par["SCvel_SSB_ini"]
    return dict(DYN_TRA_out)
