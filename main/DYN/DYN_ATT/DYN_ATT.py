# Level 2 Module DYN_ATT
# Simulates the propagation of the spacecraft attitude for the simulation

import copy
import numpy as np

from Utils import quaternions
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
            "BOFq_SSB" : par["BOFq_SSB_ini"],
            "SCw_BOF"  : par["SCw_BOF_ini"]
        }
        super().__init__("DYN_ATT", par)

    # Initalization
    def initialize(self, parent_states, DYN_states):
        # Let DYN_ATT be integrated
        self.is_dynamic = True

        self.last_input_attitude = self.par["BOFq_SSB_ini"]
        self.state["BOFq_SSB"]   = self.par["BOFq_SSB_ini"]
        self.state["SCw_BOF"]    = self.par["SCw_BOF_ini"]
        return self.state

    # Module main function
    def update_algebraic(self, t, parent_states, DYN_states, inputs):
        inputs_BOFq_SSB = inputs["DYN"]["DYN_ATT"]["BOFq_SSB"]
        inputs_SCw_BOF  = inputs["DYN"]["DYN_ATT"]["SCw_BOF"]

        # Update by orientation
        if (inputs_BOFq_SSB != self.last_input_attitude).all():
            self.state["BOFq_SSB"]   = quaternions.qnorm(inputs_BOFq_SSB)
            self.last_input_attitude = inputs_BOFq_SSB
        # Update by angular speed
        else:
            self.state["SCw_BOF"] = inputs_SCw_BOF

        return self.state

    def derivatives(self, t, DYN_states):
        BOFq_SSB  = self.state["BOFq_SSB"]
        SCw_BOF   = self.state["SCw_BOF"]
        SCw_BOF_q = np.hstack([0.0, SCw_BOF])
        dq        = 0.5 * quaternions.qprod(BOFq_SSB, SCw_BOF_q)
        return dq

    # Return integrated variables
    def get_state(self):
        return self.state["BOFq_SSB"]

    # Update integrated variables into the state dict
    def set_state(self, vec):
        self.state["BOFq_SSB"] = quaternions.qnorm(vec)
        return self.state
