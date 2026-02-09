# Level 2 Module DYN_PSR
# Simulates the pulsar directions and timing relative to the SSB frame

import copy
import numpy as np

import psrqpy

from Utils.level2module import Level2Module
from .DYN_PSR_par import DYN_PSR_par

class DYN_PSR(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_PSR_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Dummy state
        self.state = {
            "PULSARSdir_SSB" : par["PULSARSdir_SSB_ini"]
        }
        super().__init__("DYN_PSR", par)

    # Initialization
    def initialize(self, states):
        # Query Pulsar Database
        # http://www.atnf.csiro.au/research/pulsar/psrcat/
        query = psrqpy.QueryATNF(params=["RAJD","DECJD","F0","F1","F2"], condition="f0 > 50 && TYPE(HE) && !TYPE(BINARY)")

        f   = query["F0"]    # [Hz]     Pulse frequency
        df  = query["F1"]    # [Hz/s]   First derivative of pulse frequency
        ddf = query["F2"]    # [Hz/s/s] Second derivative of pulse frequency
        ra  = query["RAJD"]  # [deg]    Rigth Ascension
        dec = query["DECJD"] # [deg]    Declination

        # Update parameters
        PULSARSdir_SSB = self.radec2dir(ra, dec)
        PULSARSdir_SSB = PULSARSdir_SSB.T
        self.state["PULSARSdir_SSB"] = PULSARSdir_SSB

        return self.state

    # Convert Right Ascension and Declination in Degrees to a direction vector
    def radec2dir(self, ra, dec):
        ra  = np.deg2rad(ra)
        dec = np.deg2rad(dec)

        x = np.cos(dec)*np.cos(ra)
        y = np.cos(dec)*np.sin(ra)
        z = np.sin(dec)

        return np.array([x, y, z])
