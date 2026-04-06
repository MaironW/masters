# Level 2 Module SEN_CMB
# Simulates the temperature measurements for CMBR of the moving spacecraft

import copy
import numpy as np

from Utils import quaternions
from Utils.constants import CONSTANTS_par
from Utils.level2module import Level2Module
from .SEN_CMB_par import SEN_CMB_par

class SEN_CMB(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(SEN_CMB_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Dummy state
        self.state = {
            "SEN_CMBoutflg" : par["SEN_CMBoutflg_ini"],
            "time_CMB"      : par["time_CMB_ini"],
            "T_dipole_mes"  : par["T_dipole_ini"],
        }
        # Last time update for quantization
        self._last_update_time = -par["dt"]
        self._last_state = copy.deepcopy(self.state)
        super().__init__("SEN_CMB", par)

    # Initialization
    def initialize(self, DYN_states, SEN_states):
        # CMBR temperature loaded from DYN_CMB
        self.par["T_monopole"] = DYN_states["DYN_CMB"]["T_monopole"]

        # Update initial state
        self.state = self.update_algebraic(0, DYN_states, SEN_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, SEN_states, inputs=None):
        # Output flag
        if inputs != None:
            SEN_CMBoutflg = inputs["SEN"]["SEN_CMB"]["CMBenableflg"]
            self.state["SEN_CMBoutflg"] = SEN_CMBoutflg

        # Make all outputs invalid if SEN_CMBoutflg is zero
        # This simulates that the CMB detector was turned OFF
        if self.state["SEN_CMBoutflg"] == 0:
            self.state["time_CMB"]     = 0
            self.state["T_dipole_mes"] = self.par["T_dipole_mes_ini"]

        # CMB output is valid
        else:
            # Get spacecraft velocity vector
            SCvel_SSB = DYN_states["DYN_TRA"]["SCvel_SSB"]

            # Get spacecraft orientation
            BOFq_SSB = DYN_states["DYN_ATT"]["BOFq_SSB"]

            # Compute angle between velocity vector and sensor direction
            CMBq_BOF       = self.par["CMBq_BOF"]
            CMBq_SSB       = quaternions.qprod(CMBq_BOF, BOFq_SSB)
            CMBdir_SSB     = quaternions.qvecprod(CMBq_SSB, [1,0,0])
            SCvel_SSB_norm = np.linalg.norm(SCvel_SSB)
            SCvel_SSB_dir  = SCvel_SSB/SCvel_SSB_norm
            cos_angle      = np.clip(SCvel_SSB_dir @ CMBdir_SSB, -1, 1)
            angle          = np.arccos(cos_angle)

            # Compute temperature dipole due to the galaxy velocity
            # Compute temperature dipole due to the spacecraft velocity within the galaxy

            # Apply noise
            # Apply time quantization to all states
            time_OBT = SEN_states["SEN_TIME"]["time_OBT"]
            if time_OBT - self._last_update_time >= self.par["dt"]:
                self.state["time_CMB"]      = time_OBT

                self._last_update_time = time_OBT
                self._last_state = copy.deepcopy(self.state)
            else:
                self.state = copy.deepcopy(self._last_state)

        return self.state
