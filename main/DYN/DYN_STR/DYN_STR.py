# Level 2 Module DYN_STR
# Simulates the propagation of the Star position relative to the SSB frame

import copy

import spiceypy

from Utils import spice
from Utils.level2module import Level2Module
from .DYN_STR_par import DYN_STR_par

class DYN_STR(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_STR_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {}
        super().__init__("DYN_STR", par)

    # Initialization
    def initialize(self, states):
        count = spiceypy.ktotal("STAR")
        print(count)
        return self.state

    # Because the stars are loaded only once, at the start of the simulation,
    # there will be no update_algebraic function