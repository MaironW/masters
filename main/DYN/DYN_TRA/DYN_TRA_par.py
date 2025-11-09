import numpy as np
from Utils.constants import CONSTANTS_par

# Parameters for Module DYN_TRA

DYN_TRA_par = {
    # "SCpos_TER_ini" : np.array([CONSTANTS_par["EARTHradius_cst"]+400,0,0]), # [km]   Initial relative position of spacecraft on TER reference frame
    "SCpos_TER_ini" : np.array([7000,0,0]), # [km]   Initial relative position of spacecraft on TER reference frame
    
    # "SCpos_ECI_ini" : np.array([CONSTANTS_par["EARTHradius_cst"]+400,0,0]), # [km]   Initial relative position of spacecraft on ECI reference frame
    "SCpos_ECI_ini" : np.array([7000,0,0]), # [km]   Initial relative position of spacecraft on ECI reference frame
    # "SCvel_ECI_ini" : np.array([0,7.67,0]), # [km/s] Initial relative velocity of spacecraft on ECI reference frame
    "SCvel_ECI_ini" : np.array([0,7.546049108166282,0]), # [km/s] Initial relative velocity of spacecraft on ECI reference frame
    
    "SCpos_MAR_ini" : np.array([1,0,0]), # [km]   Initial relative position of spacecraft on MAR reference frame
    "SCpos_SUN_ini" : np.array([1,0,0]), # [km]   Initial relative position of spacecraft on SUN reference frame
    "SCpos_SSB_ini" : np.array([1,0,0]), # [km]   Initial relative position of spacecraft on SSB reference frame
    "SCvel_SSB_ini" : np.array([1,0,0]), # [km/s] Initial relative velocity of spacecraft on SSB reference frame
}