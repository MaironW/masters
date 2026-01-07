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
            "STARSpos_BOF_mes"  : par["BODYpos_BOF_mes_ini"],
            "BOFq_SSB_mes"      : par["BOFq_SSB_mes_ini"],
            "SUNpos_BOF_mes"    : par["BODYpos_BOF_mes_ini"],
            "EARTHpos_BOF_mes"  : par["BODYpos_BOF_mes_ini"],
            "MOONpos_BOF_mes"   : par["BODYpos_BOF_mes_ini"],
            "MARSpos_BOF_mes"   : par["BODYpos_BOF_mes_ini"],
            "DEIMOSpos_BOF_mes" : par["BODYpos_BOF_mes_ini"],
            "PHOBOSpos_BOF_mes" : par["BODYpos_BOF_mes_ini"],
        }
        super().__init__("SEN_STR", par)

    # Initialization
    def initialize(self, DYN_states, SEN_states):
        # Allocate STARSpos_BOF_mes
        STARSpos_SSB  = DYN_states["DYN_STR"]["STARSpos_SSB"]
        m, n = STARSpos_SSB.shape
        STARSpos_BOF_mes = np.zeros((m, n))
        self.state["STARSpos_BOF_mes"] = STARSpos_BOF_mes

        # Initialize other variables
        self.state = self.update_algebraic(0, DYN_states, SEN_states)
        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, SEN_states, inputs=None):
        # Time
        time_SIM = DYN_states["DYN_TIME"]["time_SIM"]
        self.update_time_STR(time_SIM)

        # Measured satellite attitude
        BOFq_SSB = DYN_states["DYN_ATT"]["BOFq_SSB"]
        self.update_attitude(BOFq_SSB)

        # Measured positions relative to BOF
        SUNpos_SSB    = DYN_states["DYN_SUN"]["SUNpos_SSB"]
        EARTHpos_SSB  = DYN_states["DYN_EARTH"]["EARTHpos_SSB"]
        MOONpos_SSB   = DYN_states["DYN_EARTH"]["MOONpos_SSB"]
        MARSpos_SSB   = DYN_states["DYN_MARS"]["MARSpos_SSB"]
        DEIMOSpos_SSB = DYN_states["DYN_MARS"]["DEIMOSpos_SSB"]
        PHOBOSpos_SSB = DYN_states["DYN_MARS"]["PHOBOSpos_SSB"]
        STARSpos_SSB  = DYN_states["DYN_STR"]["STARSpos_SSB"]
        self.state["SUNpos_BOF_mes"]    = self.SSBpos_BOF_normalized(SUNpos_SSB)
        self.state["EARTHpos_BOF_mes"]  = self.SSBpos_BOF_normalized(EARTHpos_SSB)
        self.state["MOONpos_BOF_mes"]   = self.SSBpos_BOF_normalized(MOONpos_SSB)
        self.state["MARSpos_BOF_mes"]   = self.SSBpos_BOF_normalized(MARSpos_SSB)
        self.state["DEIMOSpos_BOF_mes"] = self.SSBpos_BOF_normalized(DEIMOSpos_SSB)
        self.state["PHOBOSpos_BOF_mes"] = self.SSBpos_BOF_normalized(PHOBOSpos_SSB)
        self.state["STARSpos_BOF_mes"]  = self.SSBpos_BOF_normalized(STARSpos_SSB)


        return self.state

    # Update the STR time reference based on SIM time
    def update_time_STR(self, time_SIM):
        dt = self.par["dt"]
        time_STR = np.floor(time_SIM/dt) * dt
        self.state["time_STR"] = time_STR

    # Update the spacecraft attitude measurement relative to the SSB frame
    # TODO: Add noise over attitude
    # TODO: Quanticize
    def update_attitude(self, BOFq_SSB):
        BOFq_SSB_mes = BOFq_SSB
        self.state["BOFq_SSB_mes"] = BOFq_SSB_mes

    # Convert a position expressed in SSB to a measured normalized direction in BOF frame
    # TODO: Quanticize
    def SSBpos_BOF_normalized(self, pos_SSB):
        BOFq_SSB_mes = self.state["BOFq_SSB_mes"]
        pos_BOF_mes  = quaternions.qvecrot(pos_SSB, BOFq_SSB_mes)
        pos_BOF_mes  = pos_BOF_mes/np.linalg.norm(pos_BOF_mes)
        return pos_BOF_mes
