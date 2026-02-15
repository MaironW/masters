# Level 2 Module SEN_PSR
# Simulates the Time-of-Arrival measurements for pulsars by an X-ray detector
# Applies noise based on SNR

import copy
import numpy as np

from Utils.level2module import Level2Module
from .SEN_PSR_par import SEN_PSR_par

class SEN_PSR(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(SEN_PSR_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Dummy state
        self.state = {
            "SCdt_SSB_mes" : par["SCdt_SSB_mes_ini"]
        }
        super().__init__("SEN_PSR", par)

    # Initialization
    def initialize(self, states):

        # Update initial state
        self.state = self.update_algebraic(0, states)

        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, inputs=None):
        # Output flag
        if inputs != None:
            PSRoutflg = inputs["SEN"]["SEN_PSR"]["PSRenableflg"]
            self.state["PSRoutflg"] = PSRoutflg

        # Make all outputs invalid if PSRoutflg is zero
        # This simulates that the PSR detector was turned OFF

        # PSR output is valid
            # Satellite attitude
            # X-ray detector attitude
            # Pulsar directions relative to the SSB
            # Pulsar directions relative to the Spacecraft
            # Rotate to X-ray detector frame
            # Filter out objects outside the FOV

            # Compute detector noise
            # Apply noise to true SCdt_SSB for visible pulsars

            # Apply time quantization to all states (maybe not needed for PSR)

        return self.state
