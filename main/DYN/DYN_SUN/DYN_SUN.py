# Level 2 Module DYN_SUN
# Simulates the propagation of the Sun position relative to the SSB frame

import copy
import numpy as np

from Utils import spice
from Utils.level2module import Level2Module
from .DYN_SUN_par import DYN_SUN_par

class DYN_SUN(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_SUN_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "SUNpos_SSB" : par["SUNpos_SSB_ini"],
            "SUNvel_SSB" : par["SUNvel_SSB_ini"]
        }
        super().__init__("DYN_SUN", par)

    # Initialization
    def initialize(self, states):
        self.state = self.update_algebraic(0, states)
        return self.state

    # Module main function
    def update_algebraic(self, t, states, inputs=None):
        time_UTC = states["DYN_TIME"]["time_UTC"]
        SUNpos_SSB, SUNvel_SSB = spice.get_state("SUN", time_UTC)
        self.state["SUNpos_SSB"] = SUNpos_SSB
        self.state["SUNvel_SSB"] = SUNvel_SSB
        return self.state
