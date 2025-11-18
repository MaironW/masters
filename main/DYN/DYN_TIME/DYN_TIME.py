import numpy as np
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
    # Update time according to the integrator time
    DYN_out["DYN_TIME"]["time_SIM"] = DYN_TIME_par["time_UTC_ini"] + t
    DYN_out["DYN_TIME"]["time_UTC"] = DYN_TIME_par["time_UTC_ini"] + t
    return DYN_out

# Module computation of derivatives to be integrated
def derivatives(t, DYN_out):
    return np.array([])

# Return integrated variables
def get_state(DYN_TIME_out):
    return np.array([])

# Update integrated variables into the state dict
def set_state(DYN_TIME_out, vec):
    return DYN_TIME_out
