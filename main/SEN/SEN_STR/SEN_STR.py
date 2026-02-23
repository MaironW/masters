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
            "STRoutflg"         : par["STRoutflg_ini"],
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
        # Get total number of stars simulated in DYN
        STARSdir_SSB  = DYN_states["DYN_STR"]["STARSdir_SSB"]
        m, n = STARSdir_SSB.shape
        # Allocate initial stars array
        self.par["STARSid_mes_ini"]    = np.full(m, np.nan)
        self.par["STARSdir_mes_ini"]   = np.full((m,n), np.nan)
        self.state["STARSid_mes"]      = self.par["STARSid_mes_ini"]
        self.state["STARSdir_STR_mes"] = self.par["STARSdir_mes_ini"]
        self.state["STARSdir_SC_mes"]  = self.par["STARSdir_mes_ini"]

        # Initialize other variables
        self.state = self.update_algebraic(0, DYN_states, SEN_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, SEN_states, inputs=None):
        # Output flag
        if inputs != None:
            STRoutflg = inputs["SEN"]["SEN_STR"]["STRenableflg"]
            self.state["STRoutflg"] = STRoutflg

        # Make all outputs invalid if STRoutflg is zero
        # This simulates that the STR was turned OFF
        if self.state["STRoutflg"] == 0:
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

            # Stars directions relative to the Spacecraft
            STARSdir_SC = STARSdir_SSB.copy()

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
            is_nan_mask = np.isnan(STARSdir_STR).any(-1)
            STARSid_mes = np.where(is_nan_mask, np.nan, STARSid)

            # Compute noise quaternion (the same for all objects)
            noise_mean = self.par["noise_mean"]
            noise_std  = [self.par["noise_std"], self.par["noise_std"], self.par["noise_std"]]
            noise_STR  = np.random.normal(noise_mean, noise_std, size=3)
            noiseq_STR = quaternions.rotvec2q(noise_STR)

            # Apply noise on the focal plane
            SUNdir_STR_mes    = quaternions.qvecprod(noiseq_STR, SUNdir_STR)
            EARTHdir_STR_mes  = quaternions.qvecprod(noiseq_STR, EARTHdir_STR)
            MOONdir_STR_mes   = quaternions.qvecprod(noiseq_STR, MOONdir_STR)
            MARSdir_STR_mes   = quaternions.qvecprod(noiseq_STR, MARSdir_STR)
            DEIMOSdir_STR_mes = quaternions.qvecprod(noiseq_STR, DEIMOSdir_STR)
            PHOBOSdir_STR_mes = quaternions.qvecprod(noiseq_STR, PHOBOSdir_STR)
            STARSdir_STR_mes  = quaternions.qvecprod(noiseq_STR, STARSdir_STR)

            # Convert back to SC
            SSBq_STR         = quaternions.qtrans(STRq_SSB)
            SUNdir_SC_mes    = quaternions.qvecprod(SSBq_STR, SUNdir_STR_mes)
            EARTHdir_SC_mes  = quaternions.qvecprod(SSBq_STR, EARTHdir_STR_mes)
            MOONdir_SC_mes   = quaternions.qvecprod(SSBq_STR, MOONdir_STR_mes)
            MARSdir_SC_mes   = quaternions.qvecprod(SSBq_STR, MARSdir_STR_mes)
            DEIMOSdir_SC_mes = quaternions.qvecprod(SSBq_STR, DEIMOSdir_STR_mes)
            PHOBOSdir_SC_mes = quaternions.qvecprod(SSBq_STR, PHOBOSdir_STR_mes)
            STARSdir_SC_mes  = quaternions.qvecprod(SSBq_STR, STARSdir_STR_mes)

            # Compute BOFq_SSB_mes
            BOFq_STR     = quaternions.qtrans(STRq_BOF)
            SSBq_STR_mes = quaternions.qprod(SSBq_STR, noiseq_STR)
            STRq_SSB_mes = quaternions.qtrans(SSBq_STR_mes)
            BOFq_SSB_mes = quaternions.qprod(BOFq_STR, STRq_SSB_mes)

            # Apply time quantization to all states
            time_TDB = DYN_states["DYN_TIME"]["time_TDB"]
            if time_TDB - self._last_update_time >= self.par["dt"]:
                self.state["time_STR"]          = time_TDB
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

                self._last_update_time = time_TDB
                self._last_state = copy.deepcopy(self.state)
            else:
                self.state = copy.deepcopy(self._last_state)

        return self.state

    # Return the direction (unitary) vector from a position vector in any given frame
    def dir_from_pos(self, pos):
        return pos/np.linalg.norm(pos, axis=-1, keepdims=True)

    def mask_visible_objects(self, dir_STR):
        STR_boresight = np.array([0, 0, 1])
        cos_field_of_view = self.par["cos_field_of_view"]
        visible = (dir_STR @ STR_boresight) >= cos_field_of_view
        dir_STR_masked = np.full(dir_STR.shape, np.nan, dtype=float)
        dir_STR_masked[visible] = dir_STR[visible]
        return dir_STR_masked
