# Level 2 Module DYN_GRV
# Simulates the gravity model around a solar system body

import numpy as np
from .DYN_GRV_par import DYN_GRV_par
from Utils.constants import CONSTANTS_par
from Utils import quaternions

# Module output dictionary
def initialize():
    DYN_GRV_out = {
        "grvacc_TER" : DYN_GRV_par["grvacc_TER_ini"],
        "grvacc_ECI" : DYN_GRV_par["grvacc_ECI_ini"],
        "grvacc_MAR" : DYN_GRV_par["grvacc_MAR_ini"],
        "grvacc_MCI" : DYN_GRV_par["grvacc_MCI_ini"],
        "grvacc_SUN" : DYN_GRV_par["grvacc_SUN_ini"],
        "grvacc_SSB" : DYN_GRV_par["grvacc_SSB_ini"]
    }
    return DYN_GRV_out

# Module main function
def outputs(t, DYN_out):
    # Get parameters and states to make code more readable
    gravitational_cst = CONSTANTS_par["gravitational_cst"] # [km^2/kg s^2]
    SCmass_cst        = CONSTANTS_par["SCmass_cst"]        # [kg]
    EARTHmass_cst     = CONSTANTS_par["EARTHmass_cst"]     # [kg]
    SCpos_TER         = DYN_out["DYN_TRA"]["SCpos_TER"]    # [km]
    grvacc_TER        = DYN_out["DYN_GRV"]["grvacc_TER"]   # [km/s^2]

    # Compute the standard gravitational parameter around each body
    mu_EARTH = gravitational_cst*(SCmass_cst + EARTHmass_cst) # [km^3/s^2]

    # Compute the point mass acceleration (no perturbation) in the body rotative frame
    grvacc_TER = -mu_EARTH * SCpos_TER/np.linalg.norm(SCpos_TER)**3 # [km/s^2]

    # Compute the gravity acceleration in the body inertial frame
    ECIq_TER = DYN_out["DYN_EARTH"]["ECIq_TER"]
    grvacc_ECI = quaternions.qvecrot(grvacc_TER, ECIq_TER)

    DYN_out["DYN_GRV"]["grvacc_TER"] = grvacc_TER # [km/s^2]
    DYN_out["DYN_GRV"]["grvacc_ECI"] = grvacc_ECI # [km/s^2]
    return DYN_out

# Module computation of derivatives to be integrated
def derivatives(t, DYN_out):
    return np.array([])

# Return integrated variables
def get_state(DYN_GRV_out):
    return np.array([])

# Update integrated variables into the state dict
def set_state(DYN_GRV_out, vec):
    return DYN_GRV_out
