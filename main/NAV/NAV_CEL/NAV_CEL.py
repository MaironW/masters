# Level 2 Module NAV_CEL
# Provides the inputs for a Kalman Filter based on celestial measurements
# Inputs: Selected celestial bodies positions with respect to the spacecraft

import copy
import numpy as np

from Utils.level2module import Level2Module
from .NAV_CEL_par import NAV_CEL_par

class NAV_CEL(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(NAV_CEL_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "NAV_CELoutflg" : par["NAV_CELoutflg_ini"],
        }
        super().__init__("NAV_CEL", par)

    # Initialization
    def initialize(self, SEN_states, NAV_states):
        # Allocate initial arrays
        n_bodies = self.par["n_bodies"]
        self.par["z_ini"]           = np.full(n_bodies, np.nan)
        self.par["R_ini"]           = np.full((n_bodies, n_bodies), np.nan)
        self.par["STARdir_SC_ini"]  = np.full((n_bodies, 3), np.nan)
        self.par["BODYpos_SSB_ini"] = np.full((n_bodies, 3), np.nan)

        self.state["z"] = self.par["z_ini"]
        self.state["R"] = self.par["R_ini"]
        self.state["BODYsel_STARdir_SC_mes_list"] = self.par["STARdir_SC_ini"]
        self.state["BODYsel_STARdir_SC_ref_list"] = self.par["STARdir_SC_ini"]
        self.state["BODYpos_SSB_list"]            = self.par["BODYpos_SSB_ini"]
        self.state["time_valid"]                  = self.par["time_valid_ini"]

        # Update KF functions
        self.state["h"] = self.h
        self.state["H"] = self.H

        self.state = self.update_algebraic(0, SEN_states, NAV_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, SEN_states, NAV_states, inputs=None):
        SEN_STRoutflg = SEN_states["SEN_STR"]["SEN_STRoutflg"]

        # Only update outputs if SEN_STRoutflg is valid
        # Otherwise, return default values
        if SEN_STRoutflg == 0:
            self.state["NAV_CELoutflg"] = self.par["NAV_CELoutflg_ini"]
            self.state["z"]             = self.par["z_ini"]
            self.state["R"]             = self.par["R_ini"]
            self.state["time_valid"]    = self.par["time_valid_ini"]
            self.state["BODYsel_STARdir_SC_mes_list"] = self.par["STARdir_SC_ini"]
            self.state["BODYsel_STARdir_SC_ref_list"] = self.par["STARdir_SC_ini"]
            self.state["BODYpos_SSB_list"]            = self.par["BODYpos_SSB_ini"]

        # STR output is valid
        else:
            # Load sensor output
            SUNdir_SC_mes    = SEN_states["SEN_STR"]["SUNdir_SC_mes"]
            EARTHdir_SC_mes  = SEN_states["SEN_STR"]["EARTHdir_SC_mes"]
            MOONdir_SC_mes   = SEN_states["SEN_STR"]["MOONdir_SC_mes"]
            MARSdir_SC_mes   = SEN_states["SEN_STR"]["MARSdir_SC_mes"]
            PHOBOSdir_SC_mes = SEN_states["SEN_STR"]["PHOBOSdir_SC_mes"]
            DEIMOSdir_SC_mes = SEN_states["SEN_STR"]["DEIMOSdir_SC_mes"]
            STARSdir_SC_mes  = SEN_states["SEN_STR"]["STARSdir_SC_mes"]
            STARSid_mes      = SEN_states["SEN_STR"]["STARSid_mes"]
            time_STR         = SEN_states["SEN_STR"]["time_STR"]

            # Load catalog stars (assumes SSB == SC frame)
            STARSid_ref     = NAV_states["NAV_STR"]["STARSid"]
            STARSdir_SC_ref = NAV_states["NAV_STR"]["STARSdir_SSB"]

            # Filter out non-visible stars
            idx = ~np.isnan(STARSdir_SC_mes).any(axis=1)
            STARSid_mes     = STARSid_mes[idx]
            STARSdir_SC_mes = STARSdir_SC_mes[idx]

            # Check each body visibility
            SUNvisibility    = self.body_visibility(SUNdir_SC_mes)
            EARTHvisibility  = self.body_visibility(EARTHdir_SC_mes)
            MOONvisibility   = self.body_visibility(MOONdir_SC_mes)
            MARSvisibility   = self.body_visibility(MARSdir_SC_mes)
            PHOBOSvisibility = self.body_visibility(PHOBOSdir_SC_mes)
            DEIMOSvisibility = self.body_visibility(DEIMOSdir_SC_mes)

            # Get bodies ephemerides for measurement model
            SUNpos_SSB    = NAV_states["NAV_EPH"]["SUNpos_SSB"]
            EARTHpos_SSB  = NAV_states["NAV_EPH"]["EARTHpos_SSB"]
            MOONpos_SSB   = NAV_states["NAV_EPH"]["MOONpos_SSB"]
            MARSpos_SSB   = NAV_states["NAV_EPH"]["MARSpos_SSB"]
            DEIMOSpos_SSB = NAV_states["NAV_EPH"]["DEIMOSpos_SSB"]
            PHOBOSpos_SSB = NAV_states["NAV_EPH"]["PHOBOSpos_SSB"]

            # List bodies
            bodies = [
                (SUNvisibility,    SUNdir_SC_mes,    SUNpos_SSB),
                (EARTHvisibility,  EARTHdir_SC_mes,  EARTHpos_SSB),
                (MOONvisibility,   MOONdir_SC_mes,   MOONpos_SSB),
                (MARSvisibility,   MARSdir_SC_mes,   MARSpos_SSB),
                (DEIMOSvisibility, DEIMOSdir_SC_mes, DEIMOSpos_SSB),
                (PHOBOSvisibility, PHOBOSdir_SC_mes, PHOBOSpos_SSB),
            ]

            n_bodies = self.par["n_bodies"]
            z = np.full(n_bodies, np.nan)
            R = np.full(n_bodies, np.nan)
            BODYsel_STARdir_SC_mes_list = np.zeros((n_bodies, 3))
            BODYsel_STARdir_SC_ref_list = np.zeros((n_bodies, 3))
            BODYpos_SSB_list            = np.zeros((n_bodies, 3))
            sigma_angle = self.par["sigma_angle"] # [rad]
            count = 0
            for visibility, BODYdir_SC_mes, BODYpos_SSB in bodies:
                # Compute angles if body is visible
                cos_angle_mes, STARdir_SC_mes, STARid_mes, valid = self.cos_los_angle(visibility, BODYdir_SC_mes, STARSdir_SC_mes, STARSid_mes)
                # Compute the measurement model and covariance matrix
                if valid:
                    # Find matching reference star
                    ref_idx = np.where(STARSid_ref == STARid_mes)[0][0]
                    z[count] = np.arccos(cos_angle_mes)
                    R[count] = sigma_angle**2
                    BODYsel_STARdir_SC_mes_list[count] = STARdir_SC_mes
                    BODYsel_STARdir_SC_ref_list[count] = STARSdir_SC_ref[ref_idx]
                    BODYpos_SSB_list[count]            = BODYpos_SSB
                count += 1

            if np.sum(~np.isnan(z)) >= 3:
                z = z
                R = np.diag(R)
                time_valid = time_STR
                NAV_CELoutflg = 1
            else:
                z = self.par["z_ini"]
                R = np.diag(self.par["R_ini"])
                time_valid = self.par["time_valid_ini"]
                NAV_CELoutflg = 0

            # Update states
            self.state["NAV_CELoutflg"] = NAV_CELoutflg
            self.state["z"]             = z
            self.state["R"]             = R
            self.state["time_valid"]    = time_valid
            self.state["BODYsel_STARdir_SC_mes_list"] = BODYsel_STARdir_SC_mes_list
            self.state["BODYsel_STARdir_SC_ref_list"] = BODYsel_STARdir_SC_ref_list
            self.state["BODYpos_SSB_list"]            = BODYpos_SSB_list

            # Update KF functions
            self.state["h"] = self.h
            self.state["H"] = self.H

        return self.state

    # Return the celestial navigation measurement model for Kalman filtering
    def h(self, x):
        # x = [SCpos_SSB, SCvel_SSB]
        SCpos_SSB = x[0:3]

        BODYpos_SSB_list            = self.state["BODYpos_SSB_list"]
        BODYsel_STARdir_SC_ref_list = self.state["BODYsel_STARdir_SC_ref_list"]
        n_bodies                    = self.par["n_bodies"]

        h_vec = np.zeros(n_bodies)

        for i in range(n_bodies):
            # Compute the direction of the BODY from the estimated SC position, if the body is visible
            BODYpos_SSB = BODYpos_SSB_list[i]
            if np.any(BODYpos_SSB):
                # Get star direction
                # For now, do it with measurement
                # In the future, consider using ephemerides also
                STARdir_SC = BODYsel_STARdir_SC_ref_list[i]

                BODYlos_SC      = BODYpos_SSB - SCpos_SSB
                BODYlos_SC_norm = np.linalg.norm(BODYlos_SC)
                BODYdir_SC      = BODYlos_SC/BODYlos_SC_norm
                cos_angle       = np.clip(STARdir_SC @ BODYdir_SC, -1.0, 1.0)

                h_vec[i] = np.arccos(cos_angle)

        return h_vec

    # Return the Jacobian of the measurement model for EKF
    def H(self, x):
        # x = [SCpos_SSB, SCvel_SSB]
        SCpos_SSB = x[0:3]

        BODYpos_SSB_list            = self.state["BODYpos_SSB_list"]
        BODYsel_STARdir_SC_ref_list = self.state["BODYsel_STARdir_SC_ref_list"]

        n_bodies = self.par["n_bodies"]
        n_states = len(x)

        H_matrix = np.zeros((n_bodies, n_states))

        for i in range(n_bodies):
            # Compute the direction of the BODY from the estimated SC position
            BODYpos_SSB = BODYpos_SSB_list[i]
            if np.any(BODYpos_SSB):
                # Get star direction
                # For now, do it with measurement
                # In the future, consider using ephemerides also
                STARdir_SC = BODYsel_STARdir_SC_ref_list[i]

                BODYlos_SC      = BODYpos_SSB - SCpos_SSB
                BODYlos_SC_norm = np.linalg.norm(BODYlos_SC)
                BODYdir_SC      = BODYlos_SC/BODYlos_SC_norm
                cos_angle       = np.clip(STARdir_SC @ BODYdir_SC, -1.0, 1.0)

                den = np.sqrt(1 - cos_angle**2)
                # Avoid numerical blow-up on the denominator
                if den < 1e-12:
                    continue

                # Compute partial derivative
                dhdr = 1 / (BODYlos_SC_norm * den) * (STARdir_SC - cos_angle * BODYdir_SC)

                # Fill matrix
                H_matrix[i, 0:3] = dhdr

        return H_matrix

    # Check if body is visible by the star tracker
    def body_visibility(self, BODYdir_SC_mes):
        visible = ~np.isnan(BODYdir_SC_mes).any()
        return visible

    # Compute the cosine of one line-of-sight angle between the body and the best available star
    # Return the direction vectors for the selected star
    def cos_los_angle(self, BODYvisibility, BODYdir_SC_mes, STARSdir_SC_mes, STARSid_mes):
        # If body is visible
        if BODYvisibility == True:
            cos_angles = STARSdir_SC_mes @ BODYdir_SC_mes
            cos_angles = np.clip(cos_angles, -1.0, 1.0)

            # Filter-out angles smaller than los_angle_min and larger than los_angle_max
            valid_angles = (cos_angles >= self.par["cos_los_angle_min"]) & (cos_angles <= self.par["cos_los_angle_max"])
            # Check if there is still one valid angle
            if np.any(valid_angles):
                # Select valid stars
                cos_angles = cos_angles[valid_angles]
                BODYsel_STARdir_SC_mes = STARSdir_SC_mes[valid_angles]
                BODYsel_STARid_mes     = STARSid_mes[valid_angles]
                # Select the best angle
                best_idx = self.score_stars(BODYdir_SC_mes, BODYsel_STARdir_SC_mes)
                cos_angle = cos_angles[best_idx]
                BODYsel_STARdir_SC_mes = BODYsel_STARdir_SC_mes[best_idx]
                BODYsel_STARid_mes     = BODYsel_STARid_mes[best_idx]

                return cos_angle, BODYsel_STARdir_SC_mes, BODYsel_STARid_mes, True

        # If not valid
        return self.par["BODYangles_mes_ini"], self.par["BODYsel_STARdir_mes_ini"], None, False

    # Select the single most informative star relative to a body
    # Maximize angular separation from body direction
    def score_stars(self, BODYdir_SC_mes, STARSdir_SC_mes):
        best_idx = np.argmin(np.abs(STARSdir_SC_mes @ BODYdir_SC_mes))
        return best_idx
