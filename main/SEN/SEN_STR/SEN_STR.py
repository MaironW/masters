# Level 2 Module SEN_STR
# Simulates a Star Tracker, giving spacecraft attitude and relative positions of selected selestial bodies
# Inputs: The spacecraft real attitude and position of celestial bodies

import copy
import numpy as np

from Utils import quaternions
from Utils.level2module import Level2Module
from .SEN_STR_par import SEN_STR_par

class SEN_STR(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(SEN_STR_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "SEN_STRoutflg"     : par["SEN_STRoutflg_ini"],
            "time_STR"          : par["time_STR_ini"],
            "BOFq_SSB_mes"      : par["BOFq_SSB_mes_ini"],

            "SUNdir_STR_mes"    : par["BODYdir_mes_ini"],
            "EARTHdir_STR_mes"  : par["BODYdir_mes_ini"],
            "MOONdir_STR_mes"   : par["BODYdir_mes_ini"],
            "MARSdir_STR_mes"   : par["BODYdir_mes_ini"],
            "DEIMOSdir_STR_mes" : par["BODYdir_mes_ini"],
            "PHOBOSdir_STR_mes" : par["BODYdir_mes_ini"],

            "SUNdir_SC_mes"     : par["BODYdir_mes_ini"],
            "EARTHdir_SC_mes"   : par["BODYdir_mes_ini"],
            "MOONdir_SC_mes"    : par["BODYdir_mes_ini"],
            "MARSdir_SC_mes"    : par["BODYdir_mes_ini"],
            "DEIMOSdir_SC_mes"  : par["BODYdir_mes_ini"],
            "PHOBOSdir_SC_mes"  : par["BODYdir_mes_ini"],

            "STARSid_mes"       : par["STARSid_mes_ini"],
        }
        # Last time update for quantization
        self._last_update_time = -par["dt"]
        self._last_state = copy.deepcopy(self.state)
        super().__init__("SEN_STR", par)

    # Initialization
    def initialize(self, DYN_states, SEN_states):
        # Load time OBT
        time_OBT = SEN_states["SEN_TIME"]["time_OBT"]

        # Set initial value for time_OBT
        self.par["time_STR_ini"] = time_OBT
        self.state["time_STR"]   = self.par["time_STR_ini"]

        # Get total number of stars simulated in DYN
        STARSdir_SSB = DYN_states["DYN_STR"]["STARSdir_SSB"]
        m, n = STARSdir_SSB.shape
        # Allocate initial stars array
        self.par["STARSid_mes_ini"]    = np.full(m, np.nan)
        self.par["STARSdir_mes_ini"]   = np.full((m,n), np.nan)
        self.state["STARSid_mes"]      = self.par["STARSid_mes_ini"]
        self.state["STARSdir_STR_mes"] = self.par["STARSdir_mes_ini"]
        self.state["STARSdir_SC_mes"]  = self.par["STARSdir_mes_ini"]

        # Precomputed quantities
        STR2q_BOF = self.par["STR2q_STR1"]
        STR1_boresight = np.array([0, 0, 1])
        STR2_boresight = quaternions.qvecprod(STR2q_BOF, STR1_boresight) # relative to STR1
        self.par["STR1_boresight"] = STR1_boresight
        self.par["STR2_boresight"] = STR2_boresight

        # Initialize other variables
        self.state = self.update_algebraic(0, DYN_states, SEN_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, SEN_states, inputs=None):
        # Output flag
        if inputs != None:
            SEN_STRoutflg = inputs["SEN"]["SEN_STR"]["STRenableflg"]
            self.state["SEN_STRoutflg"] = SEN_STRoutflg

        # Make all outputs invalid if STRoutflg is zero
        # This simulates that the STR was turned OFF
        if self.state["SEN_STRoutflg"] == 0:
            self.state["time_STR"]          = self.par["time_STR_ini"]
            self.state["STARSid_mes"]       = self.par["STARSid_mes_ini"]
            self.state["SUNdir_STR_mes"]    = self.par["BODYdir_mes_ini"]
            self.state["EARTHdir_STR_mes"]  = self.par["BODYdir_mes_ini"]
            self.state["MOONdir_STR_mes"]   = self.par["BODYdir_mes_ini"]
            self.state["MARSdir_STR_mes"]   = self.par["BODYdir_mes_ini"]
            self.state["DEIMOSdir_STR_mes"] = self.par["BODYdir_mes_ini"]
            self.state["PHOBOSdir_STR_mes"] = self.par["BODYdir_mes_ini"]
            self.state["STARSdir_STR_mes"]  = self.par["STARSdir_mes_ini"]
            self.state["SUNdir_SC_mes"]     = self.par["BODYdir_mes_ini"]
            self.state["EARTHdir_SC_mes"]   = self.par["BODYdir_mes_ini"]
            self.state["MOONdir_SC_mes"]    = self.par["BODYdir_mes_ini"]
            self.state["MARSdir_SC_mes"]    = self.par["BODYdir_mes_ini"]
            self.state["DEIMOSdir_SC_mes"]  = self.par["BODYdir_mes_ini"]
            self.state["PHOBOSdir_SC_mes"]  = self.par["BODYdir_mes_ini"]
            self.state["STARSdir_SC_mes"]   = self.par["STARSdir_mes_ini"]
            self.state["BOFq_SSB_mes"]      = self.par["BOFq_SSB_mes_ini"]

        # STR output is valid
        else:
            # Apply time quantization to all states before starting computations
            # This is useful here because there are so many stars
            time_OBT = SEN_states["SEN_TIME"]["time_OBT"]
            if time_OBT - self._last_update_time < self.par["dt"]:
                self.state = copy.deepcopy(self._last_state)
                return self.state

            # Satellite attitude
            BOFq_SSB = DYN_states["DYN_ATT"]["BOFq_SSB"]

            # Star tracker attitude
            STRq_BOF = self.par["STRq_BOF"]
            STRq_SSB = quaternions.qprod(STRq_BOF, BOFq_SSB)

            # Positions relative to SSB
            SCpos_SSB     = DYN_states["DYN_TRA"]["SCpos_SSB"]
            SUNpos_SSB    = DYN_states["DYN_EPH"]["SUNpos_SSB"]
            EARTHpos_SSB  = DYN_states["DYN_EPH"]["EARTHpos_SSB"]
            MOONpos_SSB   = DYN_states["DYN_EPH"]["MOONpos_SSB"]
            MARSpos_SSB   = DYN_states["DYN_EPH"]["MARSpos_SSB"]
            DEIMOSpos_SSB = DYN_states["DYN_EPH"]["DEIMOSpos_SSB"]
            PHOBOSpos_SSB = DYN_states["DYN_EPH"]["PHOBOSpos_SSB"]

            # Get positions relative to the Spacecraft, expressed in SSB
            SUNpos_SC    = SUNpos_SSB    - SCpos_SSB
            EARTHpos_SC  = EARTHpos_SSB  - SCpos_SSB
            MOONpos_SC   = MOONpos_SSB   - SCpos_SSB
            MARSpos_SC   = MARSpos_SSB   - SCpos_SSB
            DEIMOSpos_SC = DEIMOSpos_SSB - SCpos_SSB
            PHOBOSpos_SC = PHOBOSpos_SSB - SCpos_SSB

            # Directions relative to SC
            SUNdir_SC    = self.dir_from_pos(SUNpos_SC)
            EARTHdir_SC  = self.dir_from_pos(EARTHpos_SC)
            MOONdir_SC   = self.dir_from_pos(MOONpos_SC)
            MARSdir_SC   = self.dir_from_pos(MARSpos_SC)
            DEIMOSdir_SC = self.dir_from_pos(DEIMOSpos_SC)
            PHOBOSdir_SC = self.dir_from_pos(PHOBOSpos_SC)

            # Stars direction relative to the SSB and IDs
            STARSdir_SSB = DYN_states["DYN_STR"]["STARSdir_SSB"]
            STARSid      = DYN_states["DYN_STR"]["STARSid"]

            # Stars directions relative to the Spacecraft (just for the ease of notation)
            STARSdir_SC = STARSdir_SSB

            # Rotate to STR frame
            SUNdir_STR    = quaternions.qvecprod(STRq_SSB, SUNdir_SC)
            EARTHdir_STR  = quaternions.qvecprod(STRq_SSB, EARTHdir_SC)
            MOONdir_STR   = quaternions.qvecprod(STRq_SSB, MOONdir_SC)
            MARSdir_STR   = quaternions.qvecprod(STRq_SSB, MARSdir_SC)
            DEIMOSdir_STR = quaternions.qvecprod(STRq_SSB, DEIMOSdir_SC)
            PHOBOSdir_STR = quaternions.qvecprod(STRq_SSB, PHOBOSdir_SC)
            STARSdir_STR  = quaternions.qvecprod(STRq_SSB, STARSdir_SC)

            # Filter out objects outside of the FOV
            SUNdir_STR    = self.mask_visible_objects(SUNdir_STR)
            EARTHdir_STR  = self.mask_visible_objects(EARTHdir_STR)
            MOONdir_STR   = self.mask_visible_objects(MOONdir_STR)
            MARSdir_STR   = self.mask_visible_objects(MARSdir_STR)
            DEIMOSdir_STR = self.mask_visible_objects(DEIMOSdir_STR)
            PHOBOSdir_STR = self.mask_visible_objects(PHOBOSdir_STR)
            STARSdir_STR  = self.mask_visible_objects(STARSdir_STR)

            # Filter out STARSid outside the FOV
            visible = ~np.isnan(STARSdir_STR).any(axis=1)
            STARSid_mes = np.where(~visible, np.nan, STARSid)

            # Compute noise quaternion (different for each objects, otherwise it would cancel out on the angles)
            noise_mean = self.par["noise_mean"]
            noise_std  = [self.par["noise_std"], self.par["noise_std"], self.par["noise_std"]]
            noiseq_STR        = quaternions.rotvec2q(np.random.normal(noise_mean, noise_std, size=3))
            noiseq_STR_SUN    = quaternions.rotvec2q(np.random.normal(noise_mean, noise_std, size=3))
            noiseq_STR_EARTH  = quaternions.rotvec2q(np.random.normal(noise_mean, noise_std, size=3))
            noiseq_STR_MOON   = quaternions.rotvec2q(np.random.normal(noise_mean, noise_std, size=3))
            noiseq_STR_MARS   = quaternions.rotvec2q(np.random.normal(noise_mean, noise_std, size=3))
            noiseq_STR_DEIMOS = quaternions.rotvec2q(np.random.normal(noise_mean, noise_std, size=3))
            noiseq_STR_PHOBOS = quaternions.rotvec2q(np.random.normal(noise_mean, noise_std, size=3))

            # Apply noise on the focal plane (different value for each body)
            SUNdir_STR_mes    = quaternions.qvecprod(noiseq_STR_SUN,    SUNdir_STR)
            EARTHdir_STR_mes  = quaternions.qvecprod(noiseq_STR_EARTH,  EARTHdir_STR)
            MOONdir_STR_mes   = quaternions.qvecprod(noiseq_STR_MOON,   MOONdir_STR)
            MARSdir_STR_mes   = quaternions.qvecprod(noiseq_STR_MARS,   MARSdir_STR)
            DEIMOSdir_STR_mes = quaternions.qvecprod(noiseq_STR_DEIMOS, DEIMOSdir_STR)
            PHOBOSdir_STR_mes = quaternions.qvecprod(noiseq_STR_PHOBOS, PHOBOSdir_STR)

            # Do not apply noise to the stars
            # This is computationally easier to do than to apply different noise to each star
            # And has no downside because each body already have their own noises
            # This part could be updated to account for star tracker misalignment
            STARSdir_STR_mes = np.copy(STARSdir_STR)

            # Convert back to SC
            SSBq_STR         = quaternions.qtrans(STRq_SSB)
            SUNdir_SC_mes    = quaternions.qvecprod(SSBq_STR, SUNdir_STR_mes)
            EARTHdir_SC_mes  = quaternions.qvecprod(SSBq_STR, EARTHdir_STR_mes)
            MOONdir_SC_mes   = quaternions.qvecprod(SSBq_STR, MOONdir_STR_mes)
            MARSdir_SC_mes   = quaternions.qvecprod(SSBq_STR, MARSdir_STR_mes)
            DEIMOSdir_SC_mes = quaternions.qvecprod(SSBq_STR, DEIMOSdir_STR_mes)
            PHOBOSdir_SC_mes = quaternions.qvecprod(SSBq_STR, PHOBOSdir_STR_mes)
            # Copying STARSdir_SC_mes from STARSdir_STR_mes is only safe because the NaN values would be in the same spots for both
            STARSdir_SC_mes = np.copy(STARSdir_STR_mes)
            STARSdir_SC_mes[visible] = quaternions.qvecprod(SSBq_STR, STARSdir_STR_mes[visible])

            # Compute BOFq_SSB_mes
            BOFq_STR     = quaternions.qtrans(STRq_BOF)
            SSBq_STR_mes = quaternions.qprod(noiseq_STR, SSBq_STR)
            STRq_SSB_mes = quaternions.qtrans(SSBq_STR_mes)
            BOFq_SSB_mes = quaternions.qprod(BOFq_STR, STRq_SSB_mes)

            # Update states
            self.state["time_STR"]          = time_OBT
            self.state["STARSid_mes"]       = STARSid_mes
            self.state["SUNdir_STR_mes"]    = SUNdir_STR_mes
            self.state["EARTHdir_STR_mes"]  = EARTHdir_STR_mes
            self.state["MOONdir_STR_mes"]   = MOONdir_STR_mes
            self.state["MARSdir_STR_mes"]   = MARSdir_STR_mes
            self.state["DEIMOSdir_STR_mes"] = DEIMOSdir_STR_mes
            self.state["PHOBOSdir_STR_mes"] = PHOBOSdir_STR_mes
            self.state["STARSdir_STR_mes"]  = STARSdir_STR_mes
            self.state["SUNdir_SC_mes"]     = SUNdir_SC_mes
            self.state["EARTHdir_SC_mes"]   = EARTHdir_SC_mes
            self.state["MOONdir_SC_mes"]    = MOONdir_SC_mes
            self.state["MARSdir_SC_mes"]    = MARSdir_SC_mes
            self.state["DEIMOSdir_SC_mes"]  = DEIMOSdir_SC_mes
            self.state["PHOBOSdir_SC_mes"]  = PHOBOSdir_SC_mes
            self.state["STARSdir_SC_mes"]   = STARSdir_SC_mes
            self.state["BOFq_SSB_mes"]      = BOFq_SSB_mes

            self._last_update_time = time_OBT
            self._last_state = copy.deepcopy(self.state)

        return self.state

    # Return the direction (unitary) vector from a position vector in any given frame
    def dir_from_pos(self, pos):
        return pos/np.linalg.norm(pos, axis=-1, keepdims=True)

    # Check if objects are inside any of the two Star Trackers FOV
    def mask_visible_objects(self, dir_STR):
        STR1_boresight = self.par["STR1_boresight"]
        STR2_boresight = self.par["STR2_boresight"]
        cos_field_of_view = self.par["cos_field_of_view"]
        cos_angle_1 = dir_STR @ STR1_boresight
        cos_angle_2 = dir_STR @ STR2_boresight
        visible = ((cos_angle_1 >= cos_field_of_view ) | (cos_angle_2 >= cos_field_of_view))
        dir_STR[~visible] = np.nan
        return dir_STR
