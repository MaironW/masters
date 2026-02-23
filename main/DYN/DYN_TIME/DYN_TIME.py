# Level 2 Module DYN_TIME
# Simulates the propagation of time for the simulation

import copy

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
        self.state = {
            "time_SIM" : par["time_SIM_ini"],
            "time_TDB" : par["time_TDB_ini"],
        }
        super().__init__("DYN_TIME", par)

    # Initialization
    def initialize(self, parent_states, DYN_states):
        self.state = {
            "time_SIM" : self.par["time_SIM_ini"],
            "time_TDB" : self.par["time_TDB_ini"],
        }
        return self.state

    # Module main function
    def update_algebraic(self, t, parent_states, DYN_states, inputs=None):
        # Update time according to the integrator time
        self.state["time_SIM"] = self.par["time_SIM_ini"] + t
        self.state["time_TDB"] = self.par["time_TDB_ini"] + t
        return self.state
