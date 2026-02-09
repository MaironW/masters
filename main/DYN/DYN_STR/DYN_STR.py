# Level 2 Module DYN_STR
# Simulates the propagation of the Star position relative to the SSB frame

import copy
import numpy as np

from skyfield.api import load, Star
from skyfield.data import hipparcos

from Utils.level2module import Level2Module
from .DYN_STR_par import DYN_STR_par

class DYN_STR(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_STR_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Dummy state
        self.state = {
            "STARSdir_SSB" : par["STARSdir_SSB_ini"]
        }
        super().__init__("DYN_STR", par)

    # Initialization
    def initialize(self, states):
        kernel_dir = "Utils/kernels/"
        with load.open(kernel_dir + "hip_main.dat") as f:
            star_df = hipparcos.load_dataframe(f)

        # Filter stars by magnitude
        magnitude_min = self.par["magnitude_min"]
        star_df = star_df[star_df["magnitude"] <= magnitude_min]

        # Sort stars by mangitude
        star_df = star_df.sort_values(by='magnitude', ascending=True)

        # Collect only the maximum number of stars allowed
        num_stars_max = self.par["num_stars_max"]
        star_df = star_df.head(num_stars_max)

        # Get star coordinates
        stars = Star.from_dataframe(star_df)
        pos   = stars._position_au
        norm  = np.linalg.norm(pos, axis=0)

        # Update parameters
        STARSdir_SSB = pos/norm
        STARSdir_SSB = STARSdir_SSB.T
        self.state["STARSdir_SSB"] = STARSdir_SSB

        return self.state

    # Because the stars are loaded only once, at the start of the simulation,
    # there will be no update_algebraic function
