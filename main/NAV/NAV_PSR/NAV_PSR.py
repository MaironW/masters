# Level 2 Module NAV_PSR
# Provides the inputs for a Kalman Filter based on Pulsar measurements
# Inputs: Selected pulsars TOAs difference from the SC to the SSB

import copy
import numpy as np

from Utils.level2module import Level2Module
from Utils.pulsar_database import PulsarDatabase
from .NAV_PSR_par import NAV_PSR_par

class NAV_PSR(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(NAV_PSR_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "NAV_PSRoutflg" : par["NAV_PSRoutflg_ini"],
            "z"             : par["z_ini"],
            "R"             : par["R_ini"],
        }
        super().__init__("NAV_PSR", par)

    # Initialization
    def initialize(self, SEN_states):
         # Load Pulsar database
        kernel_dir = "Utils/kernels/"
        pulsar_data_file = kernel_dir + "pulsar.csv"
        pulsar_data = PulsarDatabase(pulsar_data_file)

        # Pulsar parameters should not change over the simulation
        self.par["name"]           = pulsar_data.name          # Pulsar name
        self.par["n_pulsars"]      = pulsar_data.n_pulsars     # Number of pulsars
        self.par["PULSARSdir_SSB"] = pulsar_data.PULSARdir_SSB # Direction of Pulsar from SSB

        # Allocate initial values
        n_pulsars = self.par["n_pulsars"]
        self.par["z_ini"] = np.full((n_pulsars), np.nan)
        self.par["R_ini"] = np.full((n_pulsars, n_pulsars), np.nan)
        self.state["z"]     = self.par["z_ini"]
        self.state["R"]     = self.par["R_ini"]

        # Update initial state
        self.state = self.update_algebraic(0, SEN_states)

        return self.state

    # Module main function
    def update_algebraic(self, t, SEN_states, inputs=None):
        PSRoutflg = SEN_states["SEN_PSR"]["PSRoutflg"]

        # Only update outputs if PSRoutflg is valid
        # Otherwise, return default values
        if PSRoutflg == 0:
            self.state["NAV_PSRoutflg"] = self.par["NAV_PSRoutflg_ini"]
            self.state["z"]             = self.par["z_ini"]
            self.state["R"]             = self.par["R_ini"]

        # PSR output is valid
        else:
            # Measurement vector
            z = SEN_states["SEN_PSR"]["SCdt_SSB_mes"]

            # Measurement covariance
            sigma_TOA = self.par["sigma_TOA"]
            R = np.diag(sigma_TOA**2)

            # Update states
            self.state["NAV_PSRoutflg"] = PSRoutflg

            self.state["z"] = z
            self.state["R"] = R

        return self.state
