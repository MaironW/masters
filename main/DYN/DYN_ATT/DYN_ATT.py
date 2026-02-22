# Level 2 Module DYN_ATT
# Simulates the propagation of the spacecraft attitude for the simulation

import copy
import numpy as np

from Utils.level2module import Level2Module
from .DYN_ATT_par import DYN_ATT_par

class DYN_ATT(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_ATT_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "BOFq_SSB" : par["BOFq_SSB_ini"]
        }
        super().__init__("DYN_ATT", par)

    # Initalization
    def initialize(self, DYN_states):
        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, inputs):
        self.state["BOFq_SSB"] = inputs["DYN"]["DYN_ATT"]["BOFq_SSB"]
        return self.state
