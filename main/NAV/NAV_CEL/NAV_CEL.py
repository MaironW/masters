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
        }
        super().__init__("NAV_CEL", par)

    # Initialization
    def initialize(self, SEN_states, NAV_states):
        self.state = self.update_algebraic(0, SEN_states, NAV_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, SEN_states, NAV_states, inputs=None):
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
            SUNangles_mes    = self.par["BODYangles_mes_ini"] if SUNvisibility    == False else self.los_angles(SUNdir_SC_mes,    STARSdir_SC_mes)
            EARTHangles_mes  = self.par["BODYangles_mes_ini"] if EARTHvisibility  == False else self.los_angles(EARTHdir_SC_mes,  STARSdir_SC_mes)
            MOONangles_mes   = self.par["BODYangles_mes_ini"] if MOONvisibility   == False else self.los_angles(MOONdir_SC_mes,   STARSdir_SC_mes)
            MARSangles_mes   = self.par["BODYangles_mes_ini"] if MARSvisibility   == False else self.los_angles(MARSdir_SC_mes,   STARSdir_SC_mes)
            PHOBOSangles_mes = self.par["BODYangles_mes_ini"] if PHOBOSvisibility == False else self.los_angles(PHOBOSdir_SC_mes, STARSdir_SC_mes)
            DEIMOSangles_mes = self.par["BODYangles_mes_ini"] if DEIMOSvisibility == False else self.los_angles(DEIMOSdir_SC_mes, STARSdir_SC_mes)

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
            self.state["SUNangles_mes"]    = SUNangles_mes
            self.state["EARTHangles_mes"]  = EARTHangles_mes
            self.state["MOONangles_mes"]   = MOONangles_mes
            self.state["MARSangles_mes"]   = MARSangles_mes
            self.state["PHOBOSangles_mes"] = PHOBOSangles_mes
            self.state["DEIMOSangles_mes"] = DEIMOSangles_mes
            self.state["NAV_CELoutflg"]    = NAV_CELoutflg

        return self.state

    # Check if body is visible by the star tracker
    def body_visibility(self, BODYdir_SC_mes):
        visible =~ np.isnan(BODYdir_SC_mes).any()
        return visible

    # Compute three line-of-sight angles between the body and the available stars
    def los_angles(self, BODYdir_SC_mes, STARSdir_SC_mes):
        cos_angles = STARSdir_SC_mes @ BODYdir_SC_mes
        cos_angles = np.clip(cos_angles, -1.0, 1.0)
        angles     = np.arccos(cos_angles)

        # Filter-out angles smaller than los_angle_min
        valid_angles = angles >= self.par["los_angle_min"]
        angles = angles[valid_angles]

        # Select the three best angles (TODO: Define a criteria)
    
        return angles


    # Return NAV_CEL output flag as valid only if there are at least two valid bodies
    def check_valid_output(self, visibility_list):
        at_least_two_valid_bodies = sum(visibility_list) >= 2
        return at_least_two_valid_bodies
