# Level 2 Module DYN_PSR
# Simulates the pulsar directions and timing relative to the SSB frame

import copy
import numpy as np

import psrqpy

from Utils.constants import CONSTANTS_par
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
        # Load current TDB time and SC position
        time_TDB = states["DYN_TIME"]["time_TDB"]
        SCpos_SSB = states["DYN_TRA"]["SCpos_SSB"]

        # Query Pulsar Database
        # http://www.atnf.csiro.au/research/pulsar/psrcat/
        query = psrqpy.QueryATNF(params=["RAJD","DECJD","PEPOCH","F0","F1","F2"], condition=self.par["selection_criteria"])

        t0_mjd = query["PEPOCH"] # [MJD]    Epoch for frequency
        f      = query["F0"]     # [Hz]     Pulse frequency
        df     = query["F1"]     # [Hz/s]   First derivative of pulse frequency
        ddf    = query["F2"]     # [Hz/s/s] Second derivative of pulse frequency
        ra     = query["RAJD"]   # [deg]    Rigth Ascension
        dec    = query["DECJD"]  # [deg]    Declination

        # Convert epoch from MJD to TDB
        t0 = (t0_mjd - CONSTANTS_par["MJD2000epoch_relMJD_TDB_days"]) * CONSTANTS_par["day2sec_cst"]

        # Compute pulsar rotational phase at the SSB
        dt_SSB = time_TDB - t0
        phase_SSB = f*dt_SSB + 1/2*df*dt_SSB**2 + 1/6*ddf*dt_SSB**3

        # Compute pulsar direction in SSB
        PULSARSdir_SSB = self.radec2dir(ra, dec)
        PULSARSdir_SSB = PULSARSdir_SSB.T

        # Compute the time as perceived by the spacecraft
        light_speed_cst = CONSTANTS_par["light_speed_cst"] # [km/s]
        roemer_delay = (PULSARSdir_SSB @ SCpos_SSB) / light_speed_cst # [s]
        time_SC = time_TDB - roemer_delay # [s]

        # Compute pulsar rotational phase at the SC
        dt_SC = time_SC - t0
        phase_SC = f*dt_SC + 1/2*df*dt_SC**2 + 1/6*ddf*dt_SC**3

        # Update parameters
        self.state["PULSARSdir_SSB"] = PULSARSdir_SSB
        self.state["phase_SSB"]      = phase_SSB
        self.state["phase_SC"]       = phase_SC

        return self.state

    # Convert Right Ascension and Declination in Degrees to a direction vector
    def radec2dir(self, ra, dec):
        ra  = np.deg2rad(ra)
        dec = np.deg2rad(dec)

        x = np.cos(dec)*np.cos(ra)
        y = np.cos(dec)*np.sin(ra)
        z = np.sin(dec)

        return np.array([x, y, z])
