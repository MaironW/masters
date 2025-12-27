# Level 2 Module DYN_TIME
# Simulates the propagation of time for the simulation

import copy
import numpy as np

from Utils.level2module import Level2Module
from .DYN_TIME_par import DYN_TIME_par

class DYN_TIME(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_TIME_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = par
        super().__init__("DYN_TIME", par)

    # Initialization
    def initialize(self, states):
        self.state = {
            "time_SIM" : self.par["time_SIM_ini"],
            "time_UTC" : self.par["time_UTC_ini"],
        }
        return self.state

    # Module main function
    def update_algebraic(self, t, states):
        # Update time according to the integrator time
        self.state["time_SIM"] = self.par["time_UTC_ini"] + t
        self.state["time_UTC"] = self.par["time_UTC_ini"] + t
        return self.state

    # Module computation of derivatives to be integrated
    def derivatives(self, t):
        return np.array([])

    # Return integrated variables
    def get_state(self):
        return np.array([])

    # Update integrated variables into the state dict
    def set_state(self, vec):
        pass
