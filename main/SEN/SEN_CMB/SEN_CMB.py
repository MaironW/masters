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
            "SEN_CMBoutflg"     : par["SEN_CMBoutflg_ini"],
            "time_CMB"          : par["time_CMB_ini"],
            "T_dipole_CMB1_mes" : par["T_dipole_ini"],
            "T_dipole_CMB2_mes" : par["T_dipole_ini"],
            "T_dipole_CMB3_mes" : par["T_dipole_ini"],
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
            self.state["time_CMB"]          = 0
            self.state["T_dipole_CMB1_mes"] = self.par["T_dipole_ini"]
            self.state["T_dipole_CMB2_mes"] = self.par["T_dipole_ini"]
            self.state["T_dipole_CMB3_mes"] = self.par["T_dipole_ini"]

        # CMB output is valid
        else:
            # Spacecraft velocity vector within the Solar System
            SCvel_SSB = DYN_states["DYN_TRA"]["SCvel_SSB"] # [km/s]

            # Solar system velocity vector within the CMB the thermal bath expressed in the SSB frame
            SSBvel_CMB = CONSTANTS_par["SSBvel_CMB_cst"] # [km/s]

            # Add together the Solar System and spacecraft velocities
            SCvel_CMB      = SCvel_SSB + SSBvel_CMB    # [km/s]
            SCvel_CMB_norm = np.linalg.norm(SCvel_CMB) # [km/s]
            SCvel_CMB_dir  = SCvel_CMB/SCvel_CMB_norm

            # Spacecraft orientation
            BOFq_SSB = DYN_states["DYN_ATT"]["BOFq_SSB"]

            # Compute angle between velocity vector and each sensor direction
            CMB1q_BOF = self.par["CMB1q_BOF"] # Orientation of the CMB sensor with respect to the BOF frame
            CMB2q_BOF = self.par["CMB2q_BOF"] # Orientation of the CMB sensor with respect to the BOF frame
            CMB3q_BOF = self.par["CMB3q_BOF"] # Orientation of the CMB sensor with respect to the BOF frame

            CMB1q_SSB = quaternions.qprod(CMB1q_BOF, BOFq_SSB)
            CMB2q_SSB = quaternions.qprod(CMB2q_BOF, BOFq_SSB)
            CMB3q_SSB = quaternions.qprod(CMB3q_BOF, BOFq_SSB)

            CMB1dir_SSB = quaternions.qvecprod(CMB1q_SSB, [0,0,1])
            CMB2dir_SSB = quaternions.qvecprod(CMB2q_SSB, [0,0,1])
            CMB3dir_SSB = quaternions.qvecprod(CMB3q_SSB, [0,0,1])

            cos_angle1 = np.clip(SCvel_CMB_dir @ CMB1dir_SSB, -1, 1)
            cos_angle2 = np.clip(SCvel_CMB_dir @ CMB2dir_SSB, -1, 1)
            cos_angle3 = np.clip(SCvel_CMB_dir @ CMB3dir_SSB, -1, 1)

            # Compute temperature dipole due to the spacecraft velocity within the galaxy
            light_speed_cst = CONSTANTS_par["light_speed_cst"] # [km/s]
            beta            = SCvel_CMB_norm/light_speed_cst
            T_monopole      = self.par["T_monopole"] # [K]
            T_dipole_CMB1 = np.sqrt(1-beta*beta)/(1-beta*cos_angle1)*T_monopole
            T_dipole_CMB2 = np.sqrt(1-beta*beta)/(1-beta*cos_angle2)*T_monopole
            T_dipole_CMB3 = np.sqrt(1-beta*beta)/(1-beta*cos_angle3)*T_monopole

            # Apply noise
            noise_mean = self.par["noise_mean"]
            noise_std  = self.par["noise_std"]
            noise_CMB1 = np.random.normal(noise_mean, noise_std)
            noise_CMB2 = np.random.normal(noise_mean, noise_std)
            noise_CMB3 = np.random.normal(noise_mean, noise_std)

            T_dipole_CMB1_mes = T_dipole_CMB1 + noise_CMB1
            T_dipole_CMB2_mes = T_dipole_CMB2 + noise_CMB2
            T_dipole_CMB3_mes = T_dipole_CMB3 + noise_CMB3

            # Apply time quantization to all states
            time_OBT = SEN_states["SEN_TIME"]["time_OBT"]
            if time_OBT - self._last_update_time >= self.par["dt"]:
                self.state["time_CMB"]     = time_OBT
                self.state["T_dipole_CMB1_mes"] = T_dipole_CMB1_mes
                self.state["T_dipole_CMB2_mes"] = T_dipole_CMB2_mes
                self.state["T_dipole_CMB3_mes"] = T_dipole_CMB3_mes

                self._last_update_time = time_OBT
                self._last_state = copy.deepcopy(self.state)
            else:
                self.state = copy.deepcopy(self._last_state)

        return self.state
