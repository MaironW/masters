from Utils.constants import CONSTANTS_par

# Level 2 Module DYN_TIME
# Simulates the propagation of time for the simulation

from .DYN_TIME_par import DYN_TIME_par

# Module output dictionary
DYN_TIME_out = {
    "dt"       : DYN_TIME_par["dt"],
    "time_SIM" : DYN_TIME_par["time_SIM_ini"],
    "time_UTC" : DYN_TIME_par["time_UTC_ini"]
}

# Module main function
def run():
    DYN_TIME_out["time_SIM"] += DYN_TIME_par["dt"]
    DYN_TIME_out["time_UTC"] += DYN_TIME_par["dt"]
    return dict(DYN_TIME_out)
