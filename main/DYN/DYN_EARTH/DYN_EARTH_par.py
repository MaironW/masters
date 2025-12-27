import numpy as np

# Parameters for Module DYN_EARTH

DYN_EARTH_par = {
    "EARTHpos_SSB_ini" : np.array([0,0,0]), # [km] Initial position of the Earth relative to the SSB at the initial epoch
    "EARTHvel_SSB_ini" : np.array([0,0,0]), # [km] Initial velocity of the Earth relative to the SSB at the initial epoch
    "MOONpos_SSB_ini"  : np.array([0,0,0]), # [km] Initial position of the Moon relative to the SSB at the initial epoch
    "MOONvel_SSB_ini"  : np.array([0,0,0]), # [km/s] Initial velocity of the Moon relative to the SSB at the initial epoch
}