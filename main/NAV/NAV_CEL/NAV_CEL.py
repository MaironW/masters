# Level 2 Module NAV_CEL
# Provides the spacecraft position and orientation with respect to the SSB using Celestial Navigation
# Inputs: Selected celestial bodies positions with respect to the spacecraft

import copy

from Utils.level2module import Level2Module
from .NAV_CEL_par import NAV_CEL_par

class NAV_CEL(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(NAV_CEL_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "NAV_CELoutflg" : par["NAV_CELoutflg_ini"],
        }
        super().__init__("NAV_CEL", par)

    # Initialization
    def initialize(self, SEN_states, NAV_states):
        self.state = self.update_algebraic(0, SEN_states, NAV_states)
        return self.state
    
    # Module main function
    def update_algebraic(self, t, SEN_states, NAV_states, inputs=None):
        STRoutflg = SEN_states["SEN_STR"]["STRoutflg"]

        # Only update outputs if STRoutflg is valid
        if STRoutflg == 0:
            # Handle error
            self.state["NAV_CELoutflg"] = self.par["NAV_CELoutflg_ini"]
        else:
            # 
            self.state["NAV_CELoutflg"] = 1

        return self.state
