# Level 2 Module NAV_EPH
# Simulates the propagation multiple Solar System bodies position relative to the SSB frame

import copy

from Utils import spice
from Utils.level2module import Level2Module
from .NAV_EPH_par import NAV_EPH_par

class NAV_EPH(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(NAV_EPH_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "SUNpos_SSB"    : par["SUNpos_SSB_ini"],
            "SUNvel_SSB"    : par["SUNvel_SSB_ini"],
            "EARTHpos_SSB"  : par["EARTHpos_SSB_ini"],
            "EARTHvel_SSB"  : par["EARTHvel_SSB_ini"],
            "MOONpos_SSB"   : par["MOONpos_SSB_ini"],
            "MOONvel_SSB"   : par["MOONvel_SSB_ini"],
            "MARSpos_SSB"   : par["MARSpos_SSB_ini"],
            "MARSvel_SSB"   : par["MARSvel_SSB_ini"],
            "DEIMOSpos_SSB" : par["DEIMOSpos_SSB_ini"],
            "DEIMOSvel_SSB" : par["DEIMOSvel_SSB_ini"],
            "PHOBOSpos_SSB" : par["PHOBOSpos_SSB_ini"],
            "PHOBOSvel_SSB" : par["PHOBOSvel_SSB_ini"],
        }
        super().__init__("NAV_EPH", par)

    # Initialization
    def initialize(self, SEN_states):
        self.state = self.update_algebraic(0, SEN_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, SEN_states, inputs=None):
        # Use time_STR as time reference, because it is the sensor used to detect the bodies
        # Assumes time_STR = time_TDB + bias + drift + noise
        # Because NAV_CEL will not work anyway when the STR is off, this is fine
        # However, the creation of a SEN_TIME module, exclusive to handle time_OBT is a better idea
        time_STR = SEN_states["SEN_STR"]["time_STR"]

        SUNpos_SSB,    SUNvel_SSB    = spice.get_state("SUN",    time_STR)
        EARTHpos_SSB,  EARTHvel_SSB  = spice.get_state("EARTH",  time_STR)
        MOONpos_SSB,   MOONvel_SSB   = spice.get_state("MOON",   time_STR)
        MARSpos_SSB,   MARSvel_SSB   = spice.get_state("MARS",   time_STR)
        DEIMOSpos_SSB, DEIMOSvel_SSB = spice.get_state("DEIMOS", time_STR)
        PHOBOSpos_SSB, PHOBOSvel_SSB = spice.get_state("PHOBOS", time_STR)

        self.state["SUNpos_SSB"]    = SUNpos_SSB
        self.state["SUNvel_SSB"]    = SUNvel_SSB
        self.state["EARTHpos_SSB"]  = EARTHpos_SSB
        self.state["EARTHvel_SSB"]  = EARTHvel_SSB
        self.state["MOONpos_SSB"]   = MOONpos_SSB
        self.state["MOONvel_SSB"]   = MOONvel_SSB
        self.state["MARSpos_SSB"]   = MARSpos_SSB
        self.state["MARSvel_SSB"]   = MARSvel_SSB
        self.state["DEIMOSpos_SSB"] = DEIMOSpos_SSB
        self.state["DEIMOSvel_SSB"] = DEIMOSvel_SSB
        self.state["PHOBOSpos_SSB"] = PHOBOSpos_SSB
        self.state["PHOBOSvel_SSB"] = PHOBOSvel_SSB

        return self.state
