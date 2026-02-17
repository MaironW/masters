# Level 2 Module SEN_PSR
# Simulates the Time-of-Arrival measurements for pulsars by an X-ray detector
# Applies noise based on SNR

import copy
import numpy as np

from Utils.constants import CONSTANTS_par
from Utils.level2module import Level2Module
from Utils.pulsar_database import PulsarDatabase
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
            "PSRoutflg"    : par["PSRoutflg_ini"],
            "time_PSR"     : par["time_PSR_ini"],
            "SCdt_SSB_mes" : par["SCdt_SSB_mes_ini"]
        }
        # Last time update for quantization
        self._last_update_time = -par["dt"]
        self._last_state = copy.deepcopy(self.state)
        super().__init__("SEN_PSR", par)

    # Initialization
    def initialize(self, DYN_states):
        # Load Pulsar database
        kernel_dir = "Utils/kernels/"
        pulsar_data_file = kernel_dir + "pulsar.csv"
        pulsar_data = PulsarDatabase(pulsar_data_file)

        # Pulsar parameters should not change over the simulation
        self.par["name"]      = pulsar_data.name
        self.par["n_pulsars"] = pulsar_data.n_pulsars
        self.par["f"]         = pulsar_data.f
        self.par["Fx"]        = pulsar_data.Fx
        self.par["pf"]        = pulsar_data.pf
        self.par["W"]         = pulsar_data.W
        self.par["Bx"]        = pulsar_data.Bx

        # Pulsar direction in SSB already computed by the database
        self.par["PULSARSdir_SSB"] = pulsar_data.PULSARdir_SSB

        # Compute other parameters
        self.par["P"] = 1/self.par["f"] # [s] Pulse period
        self.par["d"] = self.par["W"] / self.par["P"] # [s] Pulse duty cycle

        # Allocate initial pulsars array
        n_pulsars = self.par["n_pulsars"]
        self.par["SCdt_SSB_mes_ini"] = np.full((n_pulsars), np.nan)
        self.state["SCdt_SSB_mes"] = self.par["SCdt_SSB_mes_ini"]

        # Update initial state
        self.state = self.update_algebraic(0, DYN_states)

        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, inputs=None):
        # Output flag
        if inputs != None:
            PSRoutflg = inputs["SEN"]["SEN_PSR"]["PSRenableflg"]
            self.state["PSRoutflg"] = PSRoutflg

        # Make all outputs invalid if PSRoutflg is zero
        # This simulates that the PSR detector was turned OFF
        if self.state["PSRoutflg"] == 0:
            self.state["time_PSR"]     = 0
            self.state["SCdt_SSB_mes"] = self.par["SCdt_SSB_mes_ini"]

        # PSR output is valid
        else:
            # TODO: Filter out objects outside the sensor FOV

            # Sensor properties
            A      = self.par["detector_area"] # [m^2]
            A_cm2  = A*CONSTANTS_par["m2cm_cst"]**2
            T_obs  = self.par["dt"]
            t_bias = self.par["t_bias"]

            # X-ray properties
            n_pulsars = self.par["n_pulsars"]
            Bx        = self.par["Bx"]
            Fx        = self.par["Fx"]
            pf        = self.par["pf"]
            W         = self.par["W"]
            d         = self.par["d"]

            # Compute detector noise
            Ns_pulsed    = Fx*A_cm2*T_obs*pf       # Pulsed photon counts
            Ns_nonpulsed = Fx*A_cm2*T_obs*d*(1-pf) # Nonpulsed photon counts
            Nb           = Bx*A_cm2*T_obs*d        # Background photon counts
            SNR = Ns_pulsed / np.sqrt(Nb + Ns_nonpulsed + Ns_pulsed) # Signal to Noise Ratio
            sigma_TOA = 0.5*W / SNR
            noise = np.random.randn(n_pulsars) * sigma_TOA

            # Apply noise to true SCdt_SSB for visible pulsars
            SCdt_SSB_mes = DYN_states["DYN_PSR"]["SCdt_SSB"] + t_bias + noise

            # Apply time quantization to all states (maybe not needed for PSR)
            time_SIM = DYN_states["DYN_TIME"]["time_SIM"]
            if time_SIM - self._last_update_time >= self.par["dt"]:
                self.state["time_PSR"]     = time_SIM
                self.state["SCdt_SSB_mes"] = SCdt_SSB_mes

                self._last_update_time = time_SIM
                self._last_state = copy.deepcopy(self.state)
            else:
                self.state = copy.deepcopy(self._last_state)

        return self.state
