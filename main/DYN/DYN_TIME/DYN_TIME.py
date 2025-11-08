from Utils.constants import CONSTANTS_par

# Level 2 Module DYN_TIME
# Simulates the propagation of time for the simulation

from .DYN_TIME_par import DYN_TIME_par

# Module output dictionary
def initialize():
    DYN_TIME_out = {
        "dt"       : DYN_TIME_par["dt"],
        "time_SIM" : DYN_TIME_par["time_SIM_ini"],
        "time_UTC" : DYN_TIME_par["time_UTC_ini"]
    }
    return DYN_TIME_out

# Module main function
def outputs(t, DYN_out):
    DYN_out["DYN_TIME"]["time_SIM"] += DYN_TIME_par["dt"]
    DYN_out["DYN_TIME"]["time_UTC"] += DYN_TIME_par["dt"]
    return DYN_out
