# Level 2 Module DYN_STR
# Simulates the propagation of the Star position relative to the SSB frame

import copy
import numpy as np

from skyfield.api import load, Star
from skyfield.data import hipparcos

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
        self.state = {
            "STARSpos_SSB" : par["STARSpos_SSB_ini"]
        }
        super().__init__("DYN_STR", par)

    # Initialization
    def initialize(self, states):
        kernel_dir = "Utils/kernels/"
        with load.open(kernel_dir + "hip_main.dat") as f:
            star_df = hipparcos.load_dataframe(f)

        # Filter stars by magnitude
        star_magnitude_min = self.par["star_magnitude_min"]
        star_df = star_df[star_df["magnitude"] <= star_magnitude_min]

        # Get star coordinates
        stars = Star.from_dataframe(star_df)
        pos   = stars._position_au
        norm  = np.linalg.norm(pos, axis=0)

        # Update state
        STARSpos_SSB = pos/norm
        self.state["STARSpos_SSB"] = STARSpos_SSB

        return self.state

    # Because the stars are loaded only once, at the start of the simulation,
    # there will be no update_algebraic function