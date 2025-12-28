# Level 2 Module DYN_GRV
# Simulates the gravity model around multiple solar system bodies

import copy
import numpy as np

from Utils.constants import CONSTANTS_par
from Utils.level2module import Level2Module
from .DYN_GRV_par import DYN_GRV_par

class DYN_GRV(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_GRV_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "grvacc_SSB" : par["grvacc_SSB_ini"],
        }
        super().__init__("DYN_GRV", par)

    # Initialization
    def initialize(self, states):
        self.state = {
            "grvacc_SSB" : DYN_GRV_par["grvacc_SSB_ini"],
        }
        return self.state

    # Module main function
    def update_algebraic(self, t, states, inputs=None):
        # Get parameters and states to make code more readable
        SCpos_SCI = states["DYN_TRA"]["SCpos_SCI"] # [km]
        SCpos_ECI = states["DYN_TRA"]["SCpos_ECI"] # [km]
        SCpos_MCI = states["DYN_TRA"]["SCpos_MCI"] # [km]

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

        self.state["grvacc_SSB"] = grvacc_SSB # [km/s^2]

        return self.state
