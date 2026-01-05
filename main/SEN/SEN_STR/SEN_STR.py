# Level 2 Module SEN_STR
# Simulates a Star Tracker, giving spacecraft attitude and relative positions of selected selestial bodies
# Inputs: The spacecraft real attitude and position of celestial bodies

import copy

from Utils.level2module import Level2Module
from .SEN_STR_par import SEN_STR_par

class SEN_STR(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(SEN_STR_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {}
        super().__init__("SEN_STR", par)

    # Initialization
    def initialize(self, states):
        self.state = {}
        return self.state

    # Module main function
    def update_algebraic(self, t, states, inputs=None):
        return self.state
