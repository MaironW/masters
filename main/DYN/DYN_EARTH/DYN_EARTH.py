# Level 2 Module DYN_EARTH
# Simulates the propagation of the Earth position relative to the SSB frame

import copy
import numpy as np

from Utils import spice
from Utils.level2module import Level2Module
from .DYN_EARTH_par import DYN_EARTH_par

class DYN_EARTH(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_EARTH_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "EARTHpos_SSB" : par["EARTHpos_SSB_ini"],
            "EARTHvel_SSB" : par["EARTHvel_SSB_ini"],
            "MOONpos_SSB"  : par["MOONpos_SSB_ini"],
            "MOONvel_SSB"  : par["MOONvel_SSB_ini"],
        }
        super().__init__("DYN_EARTH", par)

    # Initialization
    def initialize(self, states):
        self.state = self.update_algebraic(0, states)
        return self.state

    # Module main function
    def update_algebraic(self, t, states, inputs=None):
        time_TDB = states["DYN_TIME"]["time_TDB"]
        EARTHpos_SSB, EARTHvel_SSB = spice.get_state("EARTH", time_TDB)
        MOONpos_SSB, MOONvel_SSB   = spice.get_state("MOON",  time_TDB)

        self.state["EARTHpos_SSB"] = EARTHpos_SSB
        self.state["EARTHvel_SSB"] = EARTHvel_SSB
        self.state["MOONpos_SSB"]  = MOONpos_SSB
        self.state["MOONvel_SSB"]  = MOONvel_SSB

        return self.state
