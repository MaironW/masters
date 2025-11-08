# Level 2 Module DYN_TRA
# Simulates the propagation of the spacecraft orbital states for the simulation

from DYN.DYN_GRV import DYN_GRV
from .DYN_TRA_par import DYN_TRA_par
from Utils.constants import CONSTANTS_par
from Utils import quaternions
import numpy as np

# Module output dictionary
def initialize():
    DYN_TRA_out = {
        "SCpos_TER" : DYN_TRA_par["SCpos_TER_ini"],
        "SCpos_ECI" : DYN_TRA_par["SCpos_ECI_ini"],
        "SCpos_MAR" : DYN_TRA_par["SCpos_MAR_ini"],
        "SCpos_SUN" : DYN_TRA_par["SCpos_SUN_ini"],
        "SCpos_SSB" : DYN_TRA_par["SCpos_SSB_ini"],
        "SCvel_ECI" : DYN_TRA_par["SCvel_ECI_ini"],
        "SCvel_SSB" : DYN_TRA_par["SCvel_SSB_ini"]
    }
    return DYN_TRA_out

# Module main function
def outputs(t, DYN_out):
    # Get parameters and states to make code more readable
    SCpos_ECI = DYN_out["DYN_TRA"]["SCpos_ECI"]  # [km]
    SCvel_ECI = DYN_out["DYN_TRA"]["SCvel_ECI"]  # [km/s]
    TERq_ECI  = DYN_out["DYN_EARTH"]["TERq_ECI"]

    SCpos_TER = quaternions.qvecrot(SCpos_ECI, TERq_ECI)

    SCpos_SSB = SCpos_ECI + DYN_out["DYN_EARTH"]["EARTHpos_SSB"]
    SCvel_SSB = SCvel_ECI + DYN_out["DYN_EARTH"]["EARTHvel_SSB"]

    DYN_out["DYN_TRA"]["SCpos_TER"] = SCpos_TER
    DYN_out["DYN_TRA"]["SCpos_ECI"] = SCpos_ECI
    DYN_out["DYN_TRA"]["SCvel_ECI"] = SCvel_ECI
    DYN_out["DYN_TRA"]["SCpos_SSB"] = SCpos_SSB
    DYN_out["DYN_TRA"]["SCvel_SSB"] = SCvel_SSB

    return DYN_out

# Module computation of derivatives to be integrated
def derivatives(t, DYN_out):
    # Get parameters and states to make code more readable
    SCvel_ECI  = DYN_out["DYN_TRA"]["SCvel_ECI"]  # [km/s]
    grvacc_ECI = DYN_out["DYN_GRV"]["grvacc_ECI"] # [km/s^2]
    # Return derivatives
    dSCpos_ECI = SCvel_ECI
    dSCvel_ECI = grvacc_ECI
    return np.hstack([dSCpos_ECI, dSCvel_ECI])

# Return integrated variables
def get_state(DYN_TRA_out):
    return np.hstack((DYN_TRA_out["SCpos_ECI"], DYN_TRA_out["SCvel_ECI"]))

# Update integrated variables into the state dict
def set_state(DYN_TRA_out, vec):
    DYN_TRA_out["SCpos_ECI"] = vec[0:3]
    DYN_TRA_out["SCvel_ECI"] = vec[3:6]
    return DYN_TRA_out