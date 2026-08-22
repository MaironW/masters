import numpy as np
from Utils import quaternions
from Utils.constants import CONSTANTS_par

# Pre-computation
field_of_view     = 180*CONSTANTS_par["deg2rad_cst"]
cos_field_of_view = np.cos(field_of_view)

# STR orientation on BOF frame
STR1q_BOF = quaternions.rotvec2q(np.array([0,-np.pi/2, 0])) # +BOFz
STR2q_BOF = quaternions.rotvec2q(np.array([np.pi/2,0, 0]))  # +BOFy
STR2q_STR1 = quaternions.qprod(STR2q_BOF, quaternions.qtrans(STR1q_BOF))

# Parameters for Module SEN_STR
SEN_STR_par = {
    # Star Tracker Parameters
    "dt"                : 60,                # [s] Discretization of the STR time
    "STRq_BOF"          : STR1q_BOF,         # Orientation of the STR1 relative to the Spacecraft body
    "STR1q_BOF"         : STR1q_BOF,         # Orientation of the STR1 relative to the Spacecraft body
    "STR2q_BOF"         : STR2q_BOF,         # Orientation of the STR2 relative to the Spacecraft body
    "STR2q_STR1"        : STR2q_STR1,        # Orientation of the STR2 relative to the STR1
    "noise_mean"        : 0,                 # [rad] STR noise mean in [X,Y,Z] directions
    "noise_std"         : 1e-5,              # [rad] STR noise std in [X,Y,Z] directions
    "field_of_view"     : field_of_view,     # [rad] STR field of view half-angle
    "cos_field_of_view" : cos_field_of_view, # Cosine of the field of view to avoid recomputations later
    # Outputs initial values
    "SEN_STRoutflg_ini" : 0,                  # STR output flag
    "time_STR_ini"      : 0,                  # [s] Star Tracker time
    "BOFq_SSB_mes_ini"  : np.full(4, np.nan), # [rad] Spacecraft orientation measured relative to the SSB frame
    "BODYdir_mes_ini"   : np.full(3, np.nan), # Generic unitary vector for a body direction
    "STARSdir_mes_ini"  : np.nan,             # Dummy, to be filled at initialization
    "STARSid_mes_ini"   : np.array([np.nan])  # Dummy, to be filled at initialization
}
