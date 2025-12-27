# Level 2 Module DYN_GRV
# Simulates the gravity model around multiple solar system bodies

import numpy as np
from .DYN_GRV_par import DYN_GRV_par
from Utils.constants import CONSTANTS_par

# Module output dictionary
def initialize(DYN_TRA_out):
    DYN_GRV_out = {
        "grvacc_SSB" : DYN_GRV_par["grvacc_SSB_ini"],
    }
    return DYN_GRV_out

# Module main function
def outputs(t, DYN_out):
    # Get parameters and states to make code more readable
    SCpos_SCI = DYN_out["DYN_TRA"]["SCpos_SCI"] # [km]
    SCpos_ECI = DYN_out["DYN_TRA"]["SCpos_ECI"] # [km]
    SCpos_MCI = DYN_out["DYN_TRA"]["SCpos_MCI"] # [km]

    # Compute the standard gravitational parameter around each body
    mu_SUN_cst   = CONSTANTS_par["mu_SUN_cst"]   # [km^3/s^2]
    mu_EARTH_cst = CONSTANTS_par["mu_EARTH_cst"] # [km^3/s^2]
    mu_MARS_cst  = CONSTANTS_par["mu_MARS_cst"]  # [km^3/s^2]

    # Compute the point mass acceleration (no perturbation) in each body inertial frame
    grvacc_SUN_SCI   = -mu_SUN_cst   * SCpos_SCI/np.linalg.norm(SCpos_SCI)**3 # [km/s^2]
    grvacc_EARTH_ECI = -mu_EARTH_cst * SCpos_ECI/np.linalg.norm(SCpos_ECI)**3 # [km/s^2]
    grvacc_MARS_MCI  = -mu_MARS_cst  * SCpos_MCI/np.linalg.norm(SCpos_MCI)**3 # [km/s^2]

    # Compute the gravity acceleration in the SSB frame (add all inertial models together)
    # Check Vallado c1.4 - Barycentric form of the N-body problem (eq 1-38)
    grvacc_SSB = grvacc_SUN_SCI + grvacc_EARTH_ECI + grvacc_MARS_MCI

    DYN_out["DYN_GRV"]["grvacc_SSB"] = grvacc_SSB # [km/s^2]

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
