# Level 2 Module NAV_CEL
# Provides the spacecraft position and orientation with respect to the SSB using Celestial Navigation
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
            "NAV_CELoutflg"    : par["NAV_CELoutflg_ini"],

            "SUNangles_mes"    : par["BODYangles_mes_ini"],
            "EARTHangles_mes"  : par["BODYangles_mes_ini"],
            "MOONangles_mes"   : par["BODYangles_mes_ini"],
            "MARSangles_mes"   : par["BODYangles_mes_ini"],
            "DEIMOSangles_mes" : par["BODYangles_mes_ini"],
            "PHOBOSangles_mes" : par["BODYangles_mes_ini"],

            "SUNsel_STARSdir_SC_mes"    : par["BODYsel_STARSdir_mes_ini"],
            "EARTHsel_STARSdir_SC_mes"  : par["BODYsel_STARSdir_mes_ini"],
            "MOONsel_STARSdir_SC_mes"   : par["BODYsel_STARSdir_mes_ini"],
            "MARSsel_STARSdir_SC_mes"   : par["BODYsel_STARSdir_mes_ini"],
            "DEIMOSsel_STARSdir_SC_mes" : par["BODYsel_STARSdir_mes_ini"],
            "PHOBOSsel_STARSdir_SC_mes" : par["BODYsel_STARSdir_mes_ini"],
        }
        super().__init__("NAV_CEL", par)

    # Initialization
    def initialize(self, SEN_states):
        self.state = self.update_algebraic(0, SEN_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, SEN_states, inputs=None):
        STRoutflg = SEN_states["SEN_STR"]["STRoutflg"]

        # Only update outputs if STRoutflg is valid
        # Otherwise, return default values
        if STRoutflg == 0:
            self.state["NAV_CELoutflg"]    = self.par["NAV_CELoutflg_ini"]
            self.state["SUNangles_mes"]    = self.par["BODYangles_mes_ini"]
            self.state["EARTHangles_mes"]  = self.par["BODYangles_mes_ini"]
            self.state["MOONangles_mes"]   = self.par["BODYangles_mes_ini"]
            self.state["MARSangles_mes"]   = self.par["BODYangles_mes_ini"]
            self.state["DEIMOSangles_mes"] = self.par["BODYangles_mes_ini"]
            self.state["PHOBOSangles_mes"] = self.par["BODYangles_mes_ini"]

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

            # Compute angles if body is visible
            SUNangles_mes,    SUNsel_STARSdir_SC_mes    = self.los_angles(SUNvisibility,    SUNdir_SC_mes,    STARSdir_SC_mes)
            EARTHangles_mes,  EARTHsel_STARSdir_SC_mes  = self.los_angles(EARTHvisibility,  EARTHdir_SC_mes,  STARSdir_SC_mes)
            MOONangles_mes,   MOONsel_STARSdir_SC_mes   = self.los_angles(MOONvisibility,   MOONdir_SC_mes,   STARSdir_SC_mes)
            MARSangles_mes,   MARSsel_STARSdir_SC_mes   = self.los_angles(MARSvisibility,   MARSdir_SC_mes,   STARSdir_SC_mes)
            DEIMOSangles_mes, DEIMOSsel_STARSdir_SC_mes = self.los_angles(DEIMOSvisibility, DEIMOSdir_SC_mes, STARSdir_SC_mes)
            PHOBOSangles_mes, PHOBOSsel_STARSdir_SC_mes = self.los_angles(PHOBOSvisibility, PHOBOSdir_SC_mes, STARSdir_SC_mes)

            # If no angle is computed, set the output flag to invalid
            visibility_list = [
                SUNvisibility,
                EARTHvisibility,
                MOONvisibility,
                MARSvisibility,
                PHOBOSvisibility,
                DEIMOSvisibility,
            ]
            NAV_CELoutflg = self.check_valid_output(visibility_list)

            # Update states
            self.state["NAV_CELoutflg"]    = NAV_CELoutflg

            self.state["SUNangles_mes"]    = SUNangles_mes
            self.state["EARTHangles_mes"]  = EARTHangles_mes
            self.state["MOONangles_mes"]   = MOONangles_mes
            self.state["MARSangles_mes"]   = MARSangles_mes
            self.state["PHOBOSangles_mes"] = PHOBOSangles_mes
            self.state["DEIMOSangles_mes"] = DEIMOSangles_mes

            self.state["SUNsel_STARSdir_SC_mes"]    = SUNsel_STARSdir_SC_mes
            self.state["EARTHsel_STARSdir_SC_mes"]  = EARTHsel_STARSdir_SC_mes
            self.state["MOONsel_STARSdir_SC_mes"]   = MOONsel_STARSdir_SC_mes
            self.state["MARSsel_STARSdir_SC_mes"]   = MARSsel_STARSdir_SC_mes
            self.state["PHOBOSsel_STARSdir_SC_mes"] = PHOBOSsel_STARSdir_SC_mes
            self.state["DEIMOSsel_STARSdir_SC_mes"] = DEIMOSsel_STARSdir_SC_mes

        return self.state

    # Check if body is visible by the star tracker
    def body_visibility(self, BODYdir_SC_mes):
        visible = ~np.isnan(BODYdir_SC_mes).any()
        return visible

    # Compute three line-of-sight angles between the body and the available stars
    # Return the direction vectors for the three selected stars
    def los_angles(self, BODYvisibility, BODYdir_SC_mes, STARSdir_SC_mes):
        # Flag to check whether there are enough valid stars
        enough_number_of_stars = True

        # If body is visible
        if BODYvisibility == True:

            cos_angles = STARSdir_SC_mes @ BODYdir_SC_mes
            cos_angles = np.clip(cos_angles, -1.0, 1.0)
            angles     = np.arccos(cos_angles)

            # Filter-out angles smaller than los_angle_min and larger than los_angle_max
            valid_angles = (angles >= self.par["los_angle_min"]) & (angles <= self.par["los_angle_max"])
            angles = angles[valid_angles]
            BODYsel_STARSdir_SC_mes = STARSdir_SC_mes[valid_angles]

            # Check if there are still three valid angles
            if len(angles) < 3:
                enough_number_of_stars = False
            else:
                # Select the three best angles
                best_idx = self.score_stars(BODYdir_SC_mes, BODYsel_STARSdir_SC_mes)
                angles = angles[best_idx]
                BODYsel_STARSdir_SC_mes = BODYsel_STARSdir_SC_mes[best_idx]

        # If body is not visible or not enough stars in view
        if BODYvisibility == False or enough_number_of_stars == False:
            angles = self.par["BODYangles_mes_ini"]
            BODYsel_STARSdir_SC_mes = self.par["BODYsel_STARSdir_mes_ini"]

        return angles, BODYsel_STARSdir_SC_mes

    # Return NAV_CEL output flag as valid only if there are at least two valid bodies
    def check_valid_output(self, visibility_list):
        at_least_two_valid_bodies = sum(visibility_list) >= 2
        return at_least_two_valid_bodies

    # Select three out of valid stars to compute angles from
    # TODO: Score by separation angle
    def score_stars(self, BODYdir_SC_mes, STARSdir_SC_mes):
        # Pick the farthest star from the body
        i = np.argmin(np.abs(STARSdir_SC_mes @ BODYdir_SC_mes))
        star1 = STARSdir_SC_mes[i]

        # Pick the farthest star from star1
        cross_norms = np.linalg.norm(np.cross(star1, STARSdir_SC_mes), axis=1)
        cross_norms[i] = -1 # exclude star1
        j = np.argmax(cross_norms)
        star2 = STARSdir_SC_mes[j]

        # Check if two first stars are too close to each other
        # Just print a warning but no stop the code
        normal_vec  = np.cross(star1, star2)
        normal_norm = np.linalg.norm(normal_vec)
        normal_dir  = normal_vec / normal_norm

        # Pick the farthest star from plane of star1 and star2
        plane_distance = np.abs(STARSdir_SC_mes @ normal_dir)
        plane_distance[[i,j]] = -1 # exclude star1 and star2
        k = np.argmax(plane_distance)

        best_idx = [i,j,k]
        return best_idx