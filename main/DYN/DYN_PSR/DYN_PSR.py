# Level 2 Module DYN_PSR
# Simulates the pulsar directions and timing relative to the SSB frame

import copy
import numpy as np

from Utils.constants import CONSTANTS_par
from Utils.level2module import Level2Module
from Utils.pulsar_database import PulsarDatabase
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
        # Load Pulsar database
        kernel_dir = "Utils/kernels/"
        pulsar_data_file = kernel_dir + "pulsar.csv"
        pulsar_data = PulsarDatabase(pulsar_data_file)

        # Pulsar parameters should not change over the simulation
        self.par["name"]  = pulsar_data.name  #        Pulsar name
        self.par["epoch"] = pulsar_data.epoch # [MJD]  Epoch for frequency
        self.par["f"]     = pulsar_data.f     # [Hz]   Pulse frequency
        self.par["df"]    = pulsar_data.df    # [Hz/s] First derivative of pulse frequency
        self.par["ra"]    = pulsar_data.ra    # [deg]  Rigth Ascension
        self.par["dec"]   = pulsar_data.dec   # [deg]  Declination

        # Pulsar direction in SSB already computed by the database
        self.state["PULSARSdir_SSB"] = pulsar_data.PULSARdir_SSB # Direction of Pulsar from SSB

        # Update initial state
        self.state = self.update_algebraic(0, states)

        return self.state

    # Module main function
    def update_algebraic(self, t, states, inputs=None):
        # Load current TDB time and SC position
        time_TDB  = states["DYN_TIME"]["time_TDB"]
        SCpos_SSB = states["DYN_TRA"]["SCpos_SSB"]

        epoch  = self.par["epoch"] # [MJD]  Epoch for frequency
        f      = self.par["f"]     # [Hz]   Pulse frequency
        df     = self.par["df"]    # [Hz/s] First derivative of pulse frequency

        # Convert epoch from MJD to TDB
        t0 = (epoch - CONSTANTS_par["MJD2000epoch_relMJD_TDB_days"]) * CONSTANTS_par["day2sec_cst"]

        # Compute pulsar rotational phase at the SSB
        dt_SSB = time_TDB - t0
        phase_SSB = self.pulsar_phase(f, df, dt_SSB)
        phase_SSB = np.mod(phase_SSB, 1.0)

        # Compute the time as perceived by the spacecraft
        PULSARSdir_SSB = self.state["PULSARSdir_SSB"]
        light_speed_cst = CONSTANTS_par["light_speed_cst"] # [km/s]
        roemer_delay = (PULSARSdir_SSB @ SCpos_SSB) / light_speed_cst # [s]
        time_SC = time_TDB - roemer_delay # [s]

        # Compute pulsar rotational phase at the SC
        dt_SC = time_SC - t0
        phase_SC = self.pulsar_phase(f, df, dt_SC)
        phase_SC = np.mod(phase_SC, 1.0)

        # Update parameters
        self.state["phase_SSB"] = phase_SSB
        self.state["phase_SC"]  = phase_SC
        self.state["SCdt_SSB"]  = roemer_delay

        return self.state

    # Return the phase of a pulsar with:
    # f:  frequency
    # df: frequency time derivative
    # dt: time that passed since epoch t0
    def pulsar_phase(self, f, df, dt):
        return f*dt + 1/2*df*dt**2
