# Level 2 Module DYN_MARS
# Simulates the propagation of the Mars position relative to the SSB frame

import copy
import numpy as np

from Utils import spice
from Utils.level2module import Level2Module
from .DYN_MARS_par import DYN_MARS_par

class DYN_MARS(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_MARS_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "MARSpos_SSB"   : par["MARSpos_SSB_ini"],
            "MARSvel_SSB"   : par["MARSvel_SSB_ini"],
            "DEIMOSpos_SSB" : par["DEIMOSpos_SSB_ini"],
            "DEIMOSvel_SSB" : par["DEIMOSvel_SSB_ini"],
            "PHOBOSpos_SSB" : par["PHOBOSpos_SSB_ini"],
            "PHOBOSvel_SSB" : par["PHOBOSvel_SSB_ini"],
        }
        super().__init__("DYN_MARS", par)

    # Initialization
    def initialize(self, states):
        self.state = self.update_algebraic(0, states)
        return self.state

    # Module main function
    def update_algebraic(self, t, states, inputs=None):
        time_UTC = states["DYN_TIME"]["time_UTC"]
        MARSpos_SSB,   MARSvel_SSB   = spice.get_state("MARS",   time_UTC)
        DEIMOSpos_SSB, DEIMOSvel_SSB = spice.get_state("DEIMOS", time_UTC)
        PHOBOSpos_SSB, PHOBOSvel_SSB = spice.get_state("PHOBOS", time_UTC)

        self.state["MARSpos_SSB"]   = MARSpos_SSB
        self.state["MARSvel_SSB"]   = MARSvel_SSB
        self.state["DEIMOSpos_SSB"] = DEIMOSpos_SSB
        self.state["DEIMOSvel_SSB"] = DEIMOSvel_SSB
        self.state["PHOBOSpos_SSB"] = PHOBOSpos_SSB
        self.state["PHOBOSvel_SSB"] = PHOBOSvel_SSB

        return self.state
