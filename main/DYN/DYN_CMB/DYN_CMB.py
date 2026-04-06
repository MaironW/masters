# Level 2 Module DYN_CMB
# Simulates the cosmic microwave background radiation as a thermal bath for an inertial observer
# The anisotropies can by added later to the monopole

import copy

from Utils.level2module import Level2Module
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
            "T_monopole" : par["T_monopole_ini"]
        }
        super().__init__("DYN_CMB", par)

    # Initalization
    def initialize(self, parent_states, DYN_states):
        return self.state
