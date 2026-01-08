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
            "STARSdir_BOF_mes"  : par["BODYdir_BOF_mes_ini"],
            "BOFq_SSB_mes"      : par["BOFq_SSB_mes_ini"],
            "SUNdir_BOF_mes"    : par["BODYdir_BOF_mes_ini"],
            "EARTHdir_BOF_mes"  : par["BODYdir_BOF_mes_ini"],
            "MOONdir_BOF_mes"   : par["BODYdir_BOF_mes_ini"],
            "MARSdir_BOF_mes"   : par["BODYdir_BOF_mes_ini"],
            "DEIMOSdir_BOF_mes" : par["BODYdir_BOF_mes_ini"],
            "PHOBOSdir_BOF_mes" : par["BODYdir_BOF_mes_ini"],
        }
        super().__init__("SEN_STR", par)

    # Initialization
    def initialize(self, DYN_states, SEN_states):
        # Allocate STARSdir_BOF_mes
        STARSdir_SSB  = DYN_states["DYN_STR"]["STARSdir_SSB"]
        m, n = STARSdir_SSB.shape
        STARSdir_BOF_mes = np.zeros((m, n))
        self.state["STARSdir_BOF_mes"] = STARSdir_BOF_mes

        # Initialize other variables
        self.state = self.update_algebraic(0, DYN_states, SEN_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, SEN_states, inputs=None):
        # Time
        time_SIM = DYN_states["DYN_TIME"]["time_SIM"]
        self.update_time_STR(time_SIM)

        # Satellite attitude
        BOFq_SSB = DYN_states["DYN_ATT"]["BOFq_SSB"]

        # Star tracker attitude
        STRq_BOF = self.par["STRq_BOF"]
        STRq_SSB = quaternions.qprod(STRq_BOF, BOFq_SSB)

        # Positions relative to SSB
        SCpos_SSB     = DYN_states["DYN_TRA"]["SCpos_SSB"]
        SUNpos_SSB    = DYN_states["DYN_SUN"]["SUNpos_SSB"]
        EARTHpos_SSB  = DYN_states["DYN_EARTH"]["EARTHpos_SSB"]
        MOONpos_SSB   = DYN_states["DYN_EARTH"]["MOONpos_SSB"]
        MARSpos_SSB   = DYN_states["DYN_MARS"]["MARSpos_SSB"]
        DEIMOSpos_SSB = DYN_states["DYN_MARS"]["DEIMOSpos_SSB"]
        PHOBOSpos_SSB = DYN_states["DYN_MARS"]["PHOBOSpos_SSB"]

        # Get positions relative to the Spacecraft, expressed in SSB
        SUNpos_SC    = SUNpos_SSB    - SCpos_SSB
        EARTHpos_SC  = EARTHpos_SSB  - SCpos_SSB
        MOONpos_SC   = MOONpos_SSB   - SCpos_SSB
        MARSpos_SC   = MARSpos_SSB   - SCpos_SSB
        DEIMOSpos_SC = DEIMOSpos_SSB - SCpos_SSB
        PHOBOSpos_SC = PHOBOSpos_SSB - SCpos_SSB

        # Directions relative to SC, expressed in the SSB frame
        SUNdir_SC    = self.dir_from_pos(SUNpos_SC)
        EARTHdir_SC  = self.dir_from_pos(EARTHpos_SC)
        MOONdir_SC   = self.dir_from_pos(MOONpos_SC)
        MARSdir_SC   = self.dir_from_pos(MARSpos_SC)
        DEIMOSdir_SC = self.dir_from_pos(DEIMOSpos_SC)
        PHOBOSdir_SC = self.dir_from_pos(PHOBOSpos_SC)

        # Stars direction relative to the SSB
        STARSdir_SSB  = DYN_states["DYN_STR"]["STARSdir_SSB"]

        # Stars directions relative to the Spacecraft
        STARSdir_SC  = STARSdir_SSB

        # Rotate to STR frame
        SUNdir_STR    = quaternions.qvecrot(SUNdir_SC,    STRq_SSB)
        EARTHdir_STR  = quaternions.qvecrot(EARTHdir_SC,  STRq_SSB)
        MOONdir_STR   = quaternions.qvecrot(MOONdir_SC,   STRq_SSB)
        MARSdir_STR   = quaternions.qvecrot(MARSdir_SC,   STRq_SSB)
        DEIMOSdir_STR = quaternions.qvecrot(DEIMOSdir_SC, STRq_SSB)
        PHOBOSdir_STR = quaternions.qvecrot(PHOBOSdir_SC, STRq_SSB)
        STARSdir_STR  = quaternions.qvecrot(STARSdir_SC,  STRq_SSB)

        # Filter out objects outside of the FOV
        self.state["SUNdir_STR"]    = self.filter_visible_objects(SUNdir_STR)
        self.state["EARTHdir_STR"]  = self.filter_visible_objects(EARTHdir_STR)
        self.state["MOONdir_STR"]   = self.filter_visible_objects(MOONdir_STR)
        self.state["MARSdir_STR"]   = self.filter_visible_objects(MARSdir_STR)
        self.state["DEIMOSdir_STR"] = self.filter_visible_objects(DEIMOSdir_STR)
        self.state["PHOBOSdir_STR"] = self.filter_visible_objects(PHOBOSdir_STR)
        self.state["STARSdir_STR"]  = self.filter_visible_objects(STARSdir_STR)

        # Compute noise quaternion

        # Apply noise to each direction vector

        # Convert back to SSB

        # Compute BOFq_SSB_mes

        # NAV_CeleNav will receive STARSdir_SSB_mes + BODYdir_SSB_mes

        return self.state

    # Update the STR time reference based on SIM time
    def update_time_STR(self, time_SIM):
        dt = self.par["dt"]
        time_STR = np.floor(time_SIM/dt) * dt
        self.state["time_STR"] = time_STR

    # Convert a position expressed in SSB to a measured normalized direction in BOF frame
    # TODO: Quanticize
    def SCpos_to_BOFdir(self, pos_SC):
        BOFq_SSB_mes = self.state["BOFq_SSB_mes"]
        dir_BOF_mes  = quaternions.qvecrot(pos_SC, BOFq_SSB_mes)
        dir_BOF_mes  = dir_BOF_mes/np.linalg.norm(dir_BOF_mes, axis=-1, keepdims=True)
        return dir_BOF_mes

    # Return the direction (unitary) vector from a position vector in any given frame
    def dir_from_pos(self, pos):
        return pos/np.linalg.norm(pos, axis=-1, keepdims=True)

    def filter_visible_objects(self, dir_STR):
        STR_boresight = np.array([0, 0, 1])
        cos_field_of_view = self.par["cos_field_of_view"]
        visible = (dir_STR @ STR_boresight) >= cos_field_of_view
        return dir_STR[visible]
