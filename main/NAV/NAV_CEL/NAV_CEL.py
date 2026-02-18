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
    def initialize(self, SEN_states):
        # Allocate initial arrays
        n_bodies = self.par["n_bodies"]
        self.par["z_ini"] = np.full(n_bodies, np.nan)
        self.par["R_ini"] = np.full((n_bodies, n_bodies), np.nan)
        self.state["z"] = self.par["z_ini"]
        self.state["R"] = self.par["R_ini"]

        self.state = self.update_algebraic(0, SEN_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, SEN_states, inputs=None):
        STRoutflg = SEN_states["SEN_STR"]["STRoutflg"]

        # Only update outputs if STRoutflg is valid
        # Otherwise, return default values
        if STRoutflg == 0:
            self.state["NAV_CELoutflg"] = self.par["NAV_CELoutflg_ini"]
            self.state["z"] = self.par["z_ini"]
            self.state["R"] = self.par["R_ini"]

        # STR output is valid
        else:
            SUNdir_SC_mes    = SEN_states["SEN_STR"]["SUNdir_SC_mes"]
            EARTHdir_SC_mes  = SEN_states["SEN_STR"]["EARTHdir_SC_mes"]
            MOONdir_SC_mes   = SEN_states["SEN_STR"]["MOONdir_SC_mes"]
            MARSdir_SC_mes   = SEN_states["SEN_STR"]["MARSdir_SC_mes"]
            PHOBOSdir_SC_mes = SEN_states["SEN_STR"]["PHOBOSdir_SC_mes"]
            DEIMOSdir_SC_mes = SEN_states["SEN_STR"]["DEIMOSdir_SC_mes"]
            STARSdir_SC_mes  = SEN_states["SEN_STR"]["STARSdir_SC_mes"]

            # Filter out non-visible stars
            STARSdir_SC_mes = STARSdir_SC_mes[~np.isnan(STARSdir_SC_mes).any(axis=1)]

            # Check each body visibility
            SUNvisibility    = self.body_visibility(SUNdir_SC_mes)
            EARTHvisibility  = self.body_visibility(EARTHdir_SC_mes)
            MOONvisibility   = self.body_visibility(MOONdir_SC_mes)
            MARSvisibility   = self.body_visibility(MARSdir_SC_mes)
            PHOBOSvisibility = self.body_visibility(PHOBOSdir_SC_mes)
            DEIMOSvisibility = self.body_visibility(DEIMOSdir_SC_mes)

            # List bodies
            bodies = [
                (SUNvisibility,    SUNdir_SC_mes),
                (EARTHvisibility,  EARTHdir_SC_mes),
                (MOONvisibility,   MOONdir_SC_mes),
                (MARSvisibility,   MARSdir_SC_mes),
                (DEIMOSvisibility, DEIMOSdir_SC_mes),
                (PHOBOSvisibility, PHOBOSdir_SC_mes),
            ]

            n_bodies = self.par["n_bodies"]
            z = np.zeros(n_bodies)
            R = np.zeros(n_bodies)
            sigma_angle = self.par["sigma_angle"] # [rad]
            count = 0
            for visibility, BODYdir_SC_mes in bodies:
                # Compute angles if body is visible
                cos_angle_mes, STARdir_SC_mes, valid = self.cos_los_angle(visibility, BODYdir_SC_mes, STARSdir_SC_mes)
                # Compute the measurement model and covariance matrix
                if valid:
                    z[count] = cos_angle_mes
                    # Variance propagation: sigma_z^2 = (1 - cos^2(angle)) * sigma_angle^2
                    R[count] = (1 - cos_angle_mes**2) * sigma_angle**2
                    count += 1

            if count >= 2:
                z = z
                R = np.diag(R)
                NAV_CELoutflg = 1
            else:
                z = self.par["z_ini"]
                R = self.par["R_ini"]
                NAV_CELoutflg = 0

            # Update states
            self.state["NAV_CELoutflg"] = NAV_CELoutflg
            self.state["z"]             = z
            self.state["R"]             = R

        return self.state

    # Check if body is visible by the star tracker
    def body_visibility(self, BODYdir_SC_mes):
        visible = ~np.isnan(BODYdir_SC_mes).any()
        return visible

    # Compute the cosine of one line-of-sight angle between the body and the best available star
    # Return the direction vectors for the selected star
    def cos_los_angle(self, BODYvisibility, BODYdir_SC_mes, STARSdir_SC_mes):
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
                # Select the best angle
                best_idx = self.score_stars(BODYdir_SC_mes, BODYsel_STARdir_SC_mes)
                cos_angle = cos_angles[best_idx]
                BODYsel_STARdir_SC_mes = BODYsel_STARdir_SC_mes[best_idx]

                return cos_angle, BODYdir_SC_mes, True

        # If not valid
        return self.par["BODYangles_mes_ini"], self.par["BODYsel_STARdir_mes_ini"], False

    # Select the single most informative star relative to a body
    # Maximize angular separation from body direction
    def score_stars(self, BODYdir_SC_mes, STARSdir_SC_mes):
        best_idx = np.argmin(np.abs(STARSdir_SC_mes @ BODYdir_SC_mes))
        return best_idx
