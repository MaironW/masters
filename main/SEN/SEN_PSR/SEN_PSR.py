# Level 2 Module SEN_PSR
# Simulates the Time-of-Arrival measurements for pulsars by an X-ray detector
# Applies noise based on SNR

import copy
import numpy as np

from Utils.level2module import Level2Module
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
    def initialize(self, states):

        # Update initial state
        self.state = self.update_algebraic(0, states)

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
            # Satellite attitude
            # X-ray detector attitude
            # Pulsar directions relative to the SSB
            # Pulsar directions relative to the Spacecraft
            # Rotate to X-ray detector frame
            # Filter out objects outside the FOV

            # Sensor properties
            A      = self.par["detector_area"]
            T_obs  = self.par["dt"]
            t_bias = self.par["t_bias"]

            # X-ray properties
            # TODO: These should be parameters only from DYN_PSR, not states
            Bx = DYN_states["DYN_PSR"]["background_photon_flux"]
            Fx = DYN_states["DYN_PSR"]["pulsar_photon_flux"]
            pf = DYN_states["DYN_PSR"]["flux_pulsed_fraction"]
            W  = DYN_states["DYN_PSR"]["pulse_width"]
            d  = DYN_states["DYN_PSR"]["pulse_duty_cycle"]

            # Compute detector noise
            Ns_total     = Fx*A*T_obs           # Total photon counts
            Ns_pulsed    = Ns_total*pf          # Pulsed photon counts
            Ns_nonpulsed = Ns_total - Ns_pulsed # Nonpulsed photon counts
            Nb           = Bx*A*T_obs           # Background photon counts
            SNR = Ns_pulsed / np.sqrt(Nb + Ns_nonpulsed + Ns_pulsed) # Signal to Noise Ratio
            sigma_TOA = 0.5*W / SNR
            n_pulsars = len(Fx)
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
