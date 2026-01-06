# Level 2 Module SEN_STR
# Simulates a Star Tracker, giving spacecraft attitude and relative positions of selected selestial bodies
# Inputs: The spacecraft real attitude and position of celestial bodies

import copy
import numpy as np

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
        self.state = {
            "STRoutflg"         : par["STRoutflg_ini"],
            "time_STR"          : par["time_STR_ini"],
            "STARSpos_BOF_mes"  : par["BODYpos_BOF_mes_ini"],
            "BOFq_SSB_mes"      : par["BOFq_SSB_mes_ini"],
            "SUNpos_BOF_mes"    : par["BODYpos_BOF_mes_ini"],
            "EARTHpos_BOF_mes"  : par["BODYpos_BOF_mes_ini"],
            "MOONpos_BOF_mes"   : par["BODYpos_BOF_mes_ini"],
            "MARSpos_BOF_mes"   : par["BODYpos_BOF_mes_ini"],
            "DEIMOSpos_BOF_mes" : par["BODYpos_BOF_mes_ini"],
            "PHOBOSpos_BOF_mes" : par["BODYpos_BOF_mes_ini"],
        }
        super().__init__("SEN_STR", par)

    # Initialization
    def initialize(self, DYN_states, SEN_states):
        self.state = self.update_algebraic(0, DYN_states, SEN_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, SEN_states, inputs=None):
        time_SIM = DYN_states["DYN_TIME"]["time_SIM"]
        self.state["time_STR"] = self.update_time_STR(time_SIM)
        return self.state

    # Update the STR time reference based on SIM time
    def update_time_STR(self, time_SIM):
        dt = self.par["dt"]
        time_STR = np.floor(time_SIM/dt) * dt
        return time_STR
