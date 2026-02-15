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
            "PULSARSdir_SSB" : par["PULSARSdir_SSB_ini"],
            "phase_SSB"      : par["phase_SSB_ini"],
            "phase_SC"       : par["phase_SC_ini"],
            "SCdt_SSB"       : par["SCdt_SSB_ini"]
        }
        super().__init__("DYN_PSR", par)

    # Initialization
    def initialize(self, states):
        # Query Pulsar Database
        # http://www.atnf.csiro.au/research/pulsar/psrcat/
        param_list = ["NAME","RAJD","DECJD","PEPOCH","F0","F1"]
        pulsar_df = psrqpy.QueryATNF(params=param_list, condition=self.par["selection_criteria"]).table.to_pandas()

        # Convert data from astropy format to numpy array
        pulsar_data = {}
        pulsar_data["NAME"]   = np.array(pulsar_df["NAME"],   dtype=str)
        pulsar_data["RAJD"]   = np.array(pulsar_df["RAJD"],   dtype=float)
        pulsar_data["DECJD"]  = np.array(pulsar_df["DECJD"],  dtype=float)
        pulsar_data["PEPOCH"] = np.array(pulsar_df["PEPOCH"], dtype=float)
        pulsar_data["F0"]     = np.array(pulsar_df["F0"],     dtype=float)
        pulsar_data["F1"]     = np.array(pulsar_df["F1"],     dtype=float)

        # Pulsar parameters should not change over the simulation
        self.par["name"]   = pulsar_data["NAME"]   #        Pulsar name
        self.par["t0_mjd"] = pulsar_data["PEPOCH"] # [MJD]  Epoch for frequency
        self.par["f"]      = pulsar_data["F0"]     # [Hz]   Pulse frequency
        self.par["df"]     = pulsar_data["F1"]     # [Hz/s] First derivative of pulse frequency
        self.par["ra"]     = pulsar_data["RAJD"]   # [deg]  Rigth Ascension
        self.par["dec"]    = pulsar_data["DECJD"]  # [deg]  Declination

        ra  = self.par["ra"]
        dec = self.par["dec"]

        # Compute pulsar direction in SSB
        PULSARSdir_SSB = self.radec2dir(ra, dec)
        self.state["PULSARSdir_SSB"] = PULSARSdir_SSB.T

        # Update initial state
        self.state = self.update_algebraic(0, states)

        return self.state

    # Module main function
    def update_algebraic(self, t, states, inputs=None):
        # Load current TDB time and SC position
        time_TDB  = states["DYN_TIME"]["time_TDB"]
        SCpos_SSB = states["DYN_TRA"]["SCpos_SSB"]

        t0_mjd = self.par["t0_mjd"] # [MJD]  Epoch for frequency
        f      = self.par["f"]      # [Hz]   Pulse frequency
        df     = self.par["df"]     # [Hz/s] First derivative of pulse frequency

        # Convert epoch from MJD to TDB
        t0 = (t0_mjd - CONSTANTS_par["MJD2000epoch_relMJD_TDB_days"]) * CONSTANTS_par["day2sec_cst"]

        # Compute pulsar rotational phase at the SSB
        dt_SSB = time_TDB - t0
        phase_SSB = f*dt_SSB + 1/2*df*dt_SSB**2
        phase_SSB = np.mod(phase_SSB, 1.0)

        # Compute the time as perceived by the spacecraft
        PULSARSdir_SSB = self.state["PULSARSdir_SSB"]
        light_speed_cst = CONSTANTS_par["light_speed_cst"] # [km/s]
        roemer_delay = (PULSARSdir_SSB @ SCpos_SSB) / light_speed_cst # [s]
        time_SC = time_TDB - roemer_delay # [s]

        # Compute pulsar rotational phase at the SC
        dt_SC = time_SC - t0
        phase_SC = f*dt_SC + 1/2*df*dt_SC**2
        phase_SC = np.mod(phase_SC, 1.0)

        # Update parameters
        self.state["phase_SSB"] = phase_SSB
        self.state["phase_SC"]  = phase_SC
        self.state["SCdt_SSB"]  = roemer_delay

        return self.state

    # Convert Right Ascension and Declination in degrees to a direction vector
    def radec2dir(self, ra, dec):
        ra  = np.deg2rad(ra)
        dec = np.deg2rad(dec)

        x = np.cos(dec)*np.cos(ra)
        y = np.cos(dec)*np.sin(ra)
        z = np.sin(dec)

        return np.array([x, y, z])
