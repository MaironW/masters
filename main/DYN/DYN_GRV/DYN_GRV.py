# Level 2 Module DYN_GRV
# Simulates the gravity model around multiple solar system bodies

import numpy as np
from .DYN_GRV_par import DYN_GRV_par
from Utils.constants import CONSTANTS_par
from Utils import quaternions

# Module output dictionary
def initialize(DYN_TRA_out):
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
    SUNmass_cst       = CONSTANTS_par["SUNmass_cst"]       # [kg]
    EARTHmass_cst     = CONSTANTS_par["EARTHmass_cst"]     # [kg]
    MARSmass_cst      = CONSTANTS_par["MARSmass_cst"]      # [kg]
    SCpos_SUN         = DYN_out["DYN_TRA"]["SCpos_SUN"]    # [km]
    SCpos_TER         = DYN_out["DYN_TRA"]["SCpos_TER"]    # [km]
    SCpos_MAR         = DYN_out["DYN_TRA"]["SCpos_MAR"]    # [km]

    # Compute the standard gravitational parameter around each body
    mu_SUN   = gravitational_cst*(SCmass_cst + SUNmass_cst)   # [km^3/s^2]
    mu_EARTH = gravitational_cst*(SCmass_cst + EARTHmass_cst) # [km^3/s^2]
    mu_MARS  = gravitational_cst*(SCmass_cst + MARSmass_cst)  # [km^3/s^2]

    # Compute the point mass acceleration (no perturbation) in each body rotative frame
    grvacc_SUN_SUN   = -mu_SUN   * SCpos_SUN/np.linalg.norm(SCpos_SUN)**3 # [km/s^2]
    grvacc_EARTH_TER = -mu_EARTH * SCpos_TER/np.linalg.norm(SCpos_TER)**3 # [km/s^2]
    grvacc_MARS_MAR  = -mu_MARS  * SCpos_MAR/np.linalg.norm(SCpos_MAR)**3 # [km/s^2]

    # Compute the gravity acceleration in each body inertial frame
    SSBq_SUN = DYN_out["DYN_SUN"]["SSBq_SUN"]
    ECIq_TER = DYN_out["DYN_EARTH"]["ECIq_TER"]
    MCIq_MAR = DYN_out["DYN_MARS"]["MCIq_MAR"]
    grvacc_SUN_SSB   = quaternions.qvecrot(grvacc_SUN_SUN, SSBq_SUN)
    grvacc_EARTH_ECI = quaternions.qvecrot(grvacc_EARTH_TER, ECIq_TER)
    grvacc_MARS_MCI  = quaternions.qvecrot(grvacc_MARS_MAR, MCIq_MAR)

    # Compute the gravity acceleration in the SSB frame (add all inertial models together)
    grvacc_SSB = grvacc_SUN_SSB + grvacc_EARTH_ECI + grvacc_MARS_MCI
    grvacc_ECI = grvacc_SSB
    grvacc_MCI = grvacc_SSB

    # Convert back the full model to each reference system
    SUNq_SSB = DYN_out["DYN_SUN"]["SUNq_SSB"]
    TERq_ECI = DYN_out["DYN_EARTH"]["TERq_ECI"]
    MARq_MCI = DYN_out["DYN_MARS"]["MARq_MCI"]    
    grvacc_SUN = quaternions.qvecrot(grvacc_SSB, SUNq_SSB)
    grvacc_TER = quaternions.qvecrot(grvacc_ECI, TERq_ECI)
    grvacc_MAR = quaternions.qvecrot(grvacc_MCI, MARq_MCI)

    DYN_out["DYN_GRV"]["grvacc_SSB"] = grvacc_SSB # [km/s^2]
    DYN_out["DYN_GRV"]["grvacc_ECI"] = grvacc_ECI # [km/s^2]
    DYN_out["DYN_GRV"]["grvacc_MCI"] = grvacc_MCI # [km/s^2]
    DYN_out["DYN_GRV"]["grvacc_SUN"] = grvacc_SUN # [km/s^2]
    DYN_out["DYN_GRV"]["grvacc_TER"] = grvacc_TER # [km/s^2]
    DYN_out["DYN_GRV"]["grvacc_MAR"] = grvacc_MAR # [km/s^2]

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
