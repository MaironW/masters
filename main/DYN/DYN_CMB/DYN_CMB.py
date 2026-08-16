# Level 2 Module DYN_CMB
# Simulates the cosmic microwave background radiation as a thermal bath for an inertial observer
# Because the anisotropies are dependent onf the observed direction, they are added in the sensor model

import copy
import numpy as np

from Utils.level2module import Level2Module
from Utils.constants import CONSTANTS_par
from .DYN_CMB_par import DYN_CMB_par

class DYN_CMB(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_CMB_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "T_monopole"    : CONSTANTS_par["T_CMBR_cst"],
            "T_isotropic"   : CONSTANTS_par["T_CMBR_cst"],
            "T_anisotropic" : 0,
        }
        super().__init__("DYN_CMB", par)

    # Initalization
    def initialize(self, parent_states, DYN_states):
        return self.state

    # Module main function
    def update_algebraic(self, t, parent_states, DYN_states, inputs=None):
        T_anisotropic = 0
        self.state["T_anisotropic"] = T_anisotropic
        self.state["T_monopole"] = self.state["T_isotropic"] + self.state["T_anisotropic"]
        return self.state
