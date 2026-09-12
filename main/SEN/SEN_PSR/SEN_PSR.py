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
            "SEN_PSRoutflg" : par["SEN_PSRoutflg_ini"],
            "time_PSR"      : par["time_PSR_ini"],
            "OBTdt_TDB_mes" : par["OBTdt_TDB_mes_ini"],
            "PULSARSid_mes" : par["PULSARSid_mes_ini"]
        }
        # Last time update for quantization
        self._last_update_time = -par["dt"]
        self._last_state = copy.deepcopy(self.state)
        super().__init__("SEN_PSR", par)

    # Initialization
    def initialize(self, DYN_states, SEN_states):
        # Load time OBT
        time_OBT = SEN_states["SEN_TIME"]["time_OBT"]

        # Set initial value for time_OBT
        self.par["time_PSR_ini"] = time_OBT
        self.state["time_PSR"] = self.par["time_PSR_ini"]

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
        self.par["OBTdt_TDB_mes_ini"] = np.full(n_pulsars, np.nan)
        self.par["PULSARSid_mes_ini"] = np.full(n_pulsars, None)
        self.state["OBTdt_TDB_mes"]   = self.par["OBTdt_TDB_mes_ini"]
        self.state["PULSARSid_mes"]   = self.par["PULSARSid_mes_ini"]

        # Sensor properties
        A       = self.par["detector_area"] # [m^2]
        A_cm2   = A*CONSTANTS_par["m2cm_cst"]**2
        T_obs   = self.par["dt"]
        SNR_max = self.par["SNR_max"]

        # X-ray properties
        n_pulsars = self.par["n_pulsars"]
        Bx        = self.par["Bx"]
        Fx        = self.par["Fx"]
        pf        = self.par["pf"]
        W         = self.par["W"]
        d         = self.par["d"]

        # Compute detector noise for all pulsars
        Ns_pulsed    = Fx*A_cm2*T_obs*pf       # Pulsed photon counts
        Ns_nonpulsed = Fx*A_cm2*T_obs*d*(1-pf) # Nonpulsed photon counts
        Nb           = Bx*A_cm2*T_obs*d        # Background photon counts
        SNR = Ns_pulsed / np.sqrt(Nb + Ns_nonpulsed + Ns_pulsed) # Signal to Noise Ratio

        # Limit SNR
        SNR = SNR_max*SNR/(SNR_max+SNR)

        # Compute std for all pulsars
        sigma_TOA = 0.5*W / SNR
        self.par["sigma_TOA"] = sigma_TOA

        # Update initial state
        self.state = self.update_algebraic(0, DYN_states, SEN_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, SEN_states, inputs=None):
        # Output flag
        if inputs != None:
            SEN_PSRoutflg = inputs["SEN"]["SEN_PSR"]["PSRenableflg"]
            self.state["SEN_PSRoutflg"] = SEN_PSRoutflg

        # Make all outputs invalid if SEN_PSRoutflg is zero
        # This simulates that the PSR detector was turned OFF
        if self.state["SEN_PSRoutflg"] == 0:
            self.state["time_PSR"]      = self.par["time_PSR_ini"]
            self.state["OBTdt_TDB_mes"] = self.par["OBTdt_TDB_mes_ini"]
            self.state["PULSARSid_mes"] = self.par["PULSARSid_mes_ini"]

        # PSR output is valid
        else:
            # Apply time quantization to all states before starting computations
            time_OBT = SEN_states["SEN_TIME"]["time_OBT"]
            if time_OBT - self._last_update_time < self.par["dt"]:
                self.state = copy.deepcopy(self._last_state)
                return self.state

            PULSARSid_mes = self.par["name"]

            # Properties
            t_bias    = self.par["t_bias"]
            n_pulsars = self.par["n_pulsars"]
            sigma_TOA = self.par["sigma_TOA"]

            # Compute noise for all pulsars
            noise = np.random.randn(n_pulsars) * sigma_TOA

            # Apply noise to true OBTdt_TDB for visible pulsars
            OBTdt_TDB_mes = DYN_states["DYN_PSR"]["OBTdt_TDB"] + t_bias + noise

            # Update states
            self.state["time_PSR"]      = time_OBT
            self.state["OBTdt_TDB_mes"] = OBTdt_TDB_mes
            self.state["PULSARSid_mes"] = PULSARSid_mes

            self._last_update_time = time_OBT
            self._last_state = copy.deepcopy(self.state)

        return self.state
