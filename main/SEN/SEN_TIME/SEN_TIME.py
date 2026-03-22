# Level 2 Module SEN_TIME
# Simulates the spacecraft onboard time

import copy
import numpy as np

from Utils.level2module import Level2Module
from .SEN_TIME_par import SEN_TIME_par

class SEN_TIME(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(SEN_TIME_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Dummy state
        self.state = {
            "time_OBT" : par["time_OBT_ini"],
        }
        # Last time update for quantization
        self._last_update_time = -par["dt"]
        self._last_state = copy.deepcopy(self.state)
        super().__init__("SEN_TIME", par)

    # Initialization
    def initialize(self, DYN_states, SEN_states):
        # Load time TDB
        time_TDB = DYN_states["DYN_TIME"]["time_TDB"]

        # Set initial value for time_OBT equal to initial time_TDB
        self.par["time_OBT_ini"] = time_TDB
        self.state = {
            "time_OBT" : self.par["time_OBT_ini"],
        }
        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, SEN_states, inputs=None):
        # Load parameters and states
        time_TDB = DYN_states["DYN_TIME"]["time_TDB"]
        sigma_t  = self.par["sigma_t"]
        bias     = self.par["bias"]

        # Compute time noise
        noise = np.random.rand() * sigma_t

        # Apply noise and bias to true measurement
        time_OBT = time_TDB + bias + noise

        # Apply time quantization to all states
        time_SIM = DYN_states["DYN_TIME"]["time_SIM"]
        if time_SIM - self._last_update_time >= self.par["dt"]:
            self.state["time_OBT"]  = time_OBT
            self._last_update_time = time_SIM
            self._last_state = copy.deepcopy(self.state)
        else:
            self.state = copy.deepcopy(self._last_state)

        return self.state
