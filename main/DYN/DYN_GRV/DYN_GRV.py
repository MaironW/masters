# Level 2 Module DYN_GRV
# Simulates the gravity model around a solar system body

import numpy as np
from .DYN_GRV_par import DYN_GRV_par
from Utils.constants import CONSTANTS_par

# Module output dictionary
DYN_GRV_out = {
    "grvacc_TER" : DYN_GRV_par["grvacc_TER_ini"],
    "grvacc_ECI" : DYN_GRV_par["grvacc_ECI_ini"],
    "grvacc_MAR" : DYN_GRV_par["grvacc_MAR_ini"],
    "grvacc_MCI" : DYN_GRV_par["grvacc_MCI_ini"],
    "grvacc_SUN" : DYN_GRV_par["grvacc_SUN_ini"],
    "grvacc_SSB" : DYN_GRV_par["grvacc_SSB_ini"]
}

# Module main function
def run(DYN_TRA_out):
    # Get parameters to make code more readable
    gravitational_cst = CONSTANTS_par["gravitational_cst"] # [km^2/kg s^2]
    SCmass_cst        = CONSTANTS_par["SCmass_cst"]        # [kg]
    SUNmass_cst       = CONSTANTS_par["SUNmass_cst"]       # [kg]
    EARTHmass_cst     = CONSTANTS_par["EARTHmass_cst"]     # [kg]
    MARSmass_cst      = CONSTANTS_par["MARSmass_cst"]      # [kg]

    # Compute the standard gravitational parameter around each body
    mu_EARTH = gravitational_cst*(SCmass_cst + EARTHmass_cst) # [km^3/s^2]
    mu_MARS  = gravitational_cst*(SCmass_cst + MARSmass_cst)  # [km^3/s^2]
    mu_SUN   = gravitational_cst*(SCmass_cst + SUNmass_cst)   # [km^3/s^2]

    # Get the distance between the spacecraft relative to each body
    SCpos_TER = DYN_TRA_out["SCpos_TER"] # [km]
    SCpos_MAR = DYN_TRA_out["SCpos_MAR"] # [km]
    SCpos_SUN = DYN_TRA_out["SCpos_SUN"] # [km]

    # Compute the point mass acceleration (no perturbation) in the body rotative frame
    grvacc_TER = -mu_EARTH * SCpos_TER/np.linalg.norm(SCpos_TER)**3 # [km/s^2]
    grvacc_MAR = -mu_MARS  * SCpos_MAR/np.linalg.norm(SCpos_MAR)**3 # [km/s^2]
    grvacc_SUN = -mu_SUN   * SCpos_SUN/np.linalg.norm(SCpos_SUN)**3 # [km/s^2]

    # Compute the gravity acceleration in the body inertial frame
    # grvacc_ECI = ECIq_TER*grvacc_TER # TODO

    DYN_GRV_out["grvacc_TER"] = grvacc_TER # [km/s^2]
    DYN_GRV_out["grvacc_MAR"] = grvacc_MAR # [km/s^2]
    DYN_GRV_out["grvacc_SUN"] = grvacc_SUN # [km/s^2]
    return dict(DYN_GRV_out)
