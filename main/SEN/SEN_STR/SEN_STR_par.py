import numpy as np
from Utils.constants import CONSTANTS_par

# Pre-computation
field_of_view     = 60*CONSTANTS_par["deg2rad_cst"]
cos_field_of_view = np.cos(field_of_view)

# Parameters for Module SEN_STR
SEN_STR_par = {
    # Star Tracker Parameters
    "dt"                : 0.1,                 # [s] Discretization of the STR time
    "STRq_BOF"          : np.array([1,0,0,0]), # Orientation of the STR relative to the Spacecraft body
    "noise_mean"        : np.array([0,0,0]),   # [rad] STR noise mean in [X,Y,Z] directions
    "noise_std"         : np.array([0,0,0]),   # [rad] STR noise var in [X,Y,Z] directions
    "field_of_view"     : field_of_view,       # [rad] STR field of view half-angle
    "cos_field_of_view" : cos_field_of_view,   # Cosine of the field of view to avoid recomputations later
    # Outputs initial values
    "STRoutflg_ini"       : 0,                          # STR output flag
    "time_STR_ini"        : 0,                          # [s] Star Tracker time
    "BOFq_SSB_mes_ini"    : np.array([1,0,0,0]),        # [rad] Spacecraft orientation measured relative to the SSB frame
    "BODYdir_BOF_mes_ini" : np.array([None,None,None]), # Generic unitary vector for a body direction
}
