# Level 2 Module NAV_CMB
# Provides the inputs for a Kalman Filter based on CMBR temperature measurements
# Inputs: CMBR temperature measurements from different spacecraft sensors

import copy
import numpy as np

from Utils import quaternions
from Utils.constants import CONSTANTS_par
from Utils.level2module import Level2Module
from .NAV_CMB_par import NAV_CMB_par

class NAV_CMB(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(NAV_CMB_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "NAV_CMBoutflg" : par["NAV_CMBoutflg_ini"],
            "z"             : par["z_ini"],
            "R"             : par["R_ini"],
            "BOFq_SSB_ref"  : par["BOFq_SSB_ini"],
        }
        super().__init__("NAV_CMB", par)

    # Initialization
    def initialize(self, SEN_states, NAV_states):
        # Update KF functions
        self.state["h"] = self.h
        self.state["H"] = self.H
        self.state["time_valid"] = self.par["time_valid_ini"]

        # Update initial state
        self.state = self.update_algebraic(0, SEN_states, NAV_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, SEN_states, NAV_states, inputs=None):
        SEN_STRoutflg  = SEN_states["SEN_STR"]["SEN_STRoutflg"]
        SEN_CMBoutflg  = SEN_states["SEN_CMB"]["SEN_CMBoutflg"]

        # Only update outputs if SEN_CMBoutflg and SEN_STRoutflg are valid
        # Otherwise, return default values
        if SEN_CMBoutflg == 0 or SEN_STRoutflg == 0:
            self.state["NAV_CMBoutflg"] = self.par["NAV_CMBoutflg_ini"]
            self.state["z"]             = self.par["z_ini"]
            self.state["R"]             = self.par["R_ini"]
            self.state["BOFq_SSB_ref"]  = self.par["BOFq_SSB_ini"]
            self.state["time_valid"]    = self.par["time_valid_ini"]

        # CMB and STR outputs are valid
        else:
            # Load SEN_STR outputs
            BOFq_SSB_ref = SEN_states["SEN_STR"]["BOFq_SSB_mes"]
            time_CMB     = SEN_states["SEN_CMB"]["time_CMB"]

            # Load SEN_CMB outputs
            T_dipole_CMB1_mes = SEN_states["SEN_CMB"]["T_dipole_CMB1_mes"] # [K]
            T_dipole_CMB2_mes = SEN_states["SEN_CMB"]["T_dipole_CMB2_mes"] # [K]
            T_dipole_CMB3_mes = SEN_states["SEN_CMB"]["T_dipole_CMB3_mes"] # [K]

            sigma_T = np.asarray(self.par["sigma_T"])
            z = np.array([T_dipole_CMB1_mes, T_dipole_CMB2_mes, T_dipole_CMB3_mes])
            R = np.diag(sigma_T**2)
            time_valid = time_CMB
            NAV_CMBoutflg = 1

            # Update states
            self.state["NAV_CMBoutflg"] = NAV_CMBoutflg
            self.state["z"]             = z
            self.state["R"]             = R
            self.state["BOFq_SSB_ref"]  = BOFq_SSB_ref
            self.state["time_valid"]    = time_valid

            # Update KF functions
            self.state["h"] = self.h
            self.state["H"] = self.H

        return self.state

    # Return the CMBR navigation measurement model for Kalman filtering
    def h(self, x):
        # x = [SCpos_SSB, SCvel_SSB]
        SCvel_SSB = x[3:6]

        # Load parameters and states
        light_speed_cst = CONSTANTS_par["light_speed_cst"] # [km/s]

        # Solar system velocity vector within the CMB the thermal bath expressed in the SSB frame
        SSBvel_CMB = CONSTANTS_par["SSBvel_CMB_cst"] # [km/s]

        # Add together the Solar System and spacecraft velocities
        SCvel_CMB      = SCvel_SSB + SSBvel_CMB    # [km/s]
        SCvel_CMB_norm = np.linalg.norm(SCvel_CMB) # [km/s]
        SCvel_CMB_dir  = SCvel_CMB/SCvel_CMB_norm

        # Spacecraft orientation as measured by SEN_STR
        BOFq_SSB_ref = self.state["BOFq_SSB_ref"]

        # Compute angle between velocity vector and each sensor direction
        CMB1q_BOF = self.par["CMB1q_BOF"] # Orientation of the CMB sensor with respect to the BOF frame
        CMB2q_BOF = self.par["CMB2q_BOF"] # Orientation of the CMB sensor with respect to the BOF frame
        CMB3q_BOF = self.par["CMB3q_BOF"] # Orientation of the CMB sensor with respect to the BOF frame

        CMB1q_SSB = quaternions.qprod(CMB1q_BOF, BOFq_SSB_ref)
        CMB2q_SSB = quaternions.qprod(CMB2q_BOF, BOFq_SSB_ref)
        CMB3q_SSB = quaternions.qprod(CMB3q_BOF, BOFq_SSB_ref)

        CMB1dir_SSB = quaternions.qvecprod(CMB1q_SSB, [0,0,1])
        CMB2dir_SSB = quaternions.qvecprod(CMB2q_SSB, [0,0,1])
        CMB3dir_SSB = quaternions.qvecprod(CMB3q_SSB, [0,0,1])

        cos_angle1 = np.clip(SCvel_CMB_dir @ CMB1dir_SSB, -1, 1)
        cos_angle2 = np.clip(SCvel_CMB_dir @ CMB2dir_SSB, -1, 1)
        cos_angle3 = np.clip(SCvel_CMB_dir @ CMB3dir_SSB, -1, 1)

        # Compute temperature dipole due to the spacecraft velocity within the galaxy
        light_speed_cst = CONSTANTS_par["light_speed_cst"] # [km/s]
        beta            = SCvel_CMB_norm/light_speed_cst
        beta            = np.clip(beta, 0.0, 1.0 - 1e-12)
        T_monopole      = self.par["T_monopole"] # [K]
        T_dipole_CMB1 = np.sqrt(1-beta*beta)/(1-beta*cos_angle1)*T_monopole
        T_dipole_CMB2 = np.sqrt(1-beta*beta)/(1-beta*cos_angle2)*T_monopole
        T_dipole_CMB3 = np.sqrt(1-beta*beta)/(1-beta*cos_angle3)*T_monopole

        h_vec = np.array([T_dipole_CMB1, T_dipole_CMB2, T_dipole_CMB3])

        return h_vec

    # Return the Jacobian of the measurment model for EKF
    def H(self, x):
        # x = [SCpos_SSB, SCvel_SSB]
        SCvel_SSB = x[3:6]

        H_matrix = np.zeros((3, 6))

        # Load parameters and states
        light_speed_cst = CONSTANTS_par["light_speed_cst"] # [km/s]

        # Solar system velocity vector within the CMB the thermal bath expressed in the SSB frame
        SSBvel_CMB = CONSTANTS_par["SSBvel_CMB_cst"] # [km/s]

        # Add together the Solar System and spacecraft velocities
        SCvel_CMB      = SCvel_SSB + SSBvel_CMB    # [km/s]
        SCvel_CMB_norm = np.linalg.norm(SCvel_CMB) # [km/s]
        SCvel_CMB_dir  = SCvel_CMB/SCvel_CMB_norm

        # Spacecraft orientation as measured by SEN_STR
        BOFq_SSB_ref = self.state["BOFq_SSB_ref"]

        # Compute angle between velocity vector and each sensor direction
        CMB1q_BOF = self.par["CMB1q_BOF"] # Orientation of the CMB sensor with respect to the BOF frame
        CMB2q_BOF = self.par["CMB2q_BOF"] # Orientation of the CMB sensor with respect to the BOF frame
        CMB3q_BOF = self.par["CMB3q_BOF"] # Orientation of the CMB sensor with respect to the BOF frame

        CMB1q_SSB = quaternions.qprod(CMB1q_BOF, BOFq_SSB_ref)
        CMB2q_SSB = quaternions.qprod(CMB2q_BOF, BOFq_SSB_ref)
        CMB3q_SSB = quaternions.qprod(CMB3q_BOF, BOFq_SSB_ref)

        CMB1dir_SSB = quaternions.qvecprod(CMB1q_SSB, [0,0,1])
        CMB2dir_SSB = quaternions.qvecprod(CMB2q_SSB, [0,0,1])
        CMB3dir_SSB = quaternions.qvecprod(CMB3q_SSB, [0,0,1])

        cos_angle1 = np.clip(SCvel_CMB_dir @ CMB1dir_SSB, -1, 1)
        cos_angle2 = np.clip(SCvel_CMB_dir @ CMB2dir_SSB, -1, 1)
        cos_angle3 = np.clip(SCvel_CMB_dir @ CMB3dir_SSB, -1, 1)

        # Compute scalars
        SCvel_SSB_norm = np.linalg.norm(SCvel_SSB)
        beta           = SCvel_SSB_norm/light_speed_cst
        beta           = np.clip(beta, 0.0, 1.0 - 1e-12)
        T_monopole     = self.par["T_monopole"] # [K]
        aux0 = np.sqrt(1 - beta*beta)
        aux1 = 1 - beta*cos_angle1
        aux2 = 1 - beta*cos_angle2
        aux3 = 1 - beta*cos_angle3

        # Compute derivative
        dh1dv = T_monopole*(SCvel_CMB_dir * -beta/(aux0*aux1*light_speed_cst)
              + aux0/aux1*aux1 * (cos_angle1*SCvel_CMB_dir/light_speed_cst + beta * (CMB1dir_SSB - cos_angle1)/SCvel_CMB_norm))

        dh2dv = T_monopole*(SCvel_CMB_dir * -beta/(aux0*aux2*light_speed_cst)
              + aux0/aux2*aux2 * (cos_angle2*SCvel_CMB_dir/light_speed_cst + beta * (CMB2dir_SSB - cos_angle2)/SCvel_CMB_norm))

        dh3dv = T_monopole*(SCvel_CMB_dir * -beta/(aux0*aux3*light_speed_cst)
              + aux0/aux3*aux3 * (cos_angle3*SCvel_CMB_dir/light_speed_cst + beta * (CMB3dir_SSB - cos_angle3)/SCvel_CMB_norm))

        # Fill matrix
        H_matrix[0, 3:6] = dh1dv
        H_matrix[1, 3:6] = dh2dv
        H_matrix[2, 3:6] = dh3dv

        return H_matrix
