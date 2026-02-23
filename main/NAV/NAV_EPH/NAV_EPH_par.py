import numpy as np

# Parameters for Module NAV_EPH

NAV_EPH_par = {
    "SUNpos_SSB_ini"    : np.array([0,0,0]), # [km]   Initial position of the Sun relative to the SSB at the initial epoch
    "SUNvel_SSB_ini"    : np.array([0,0,0]), # [km/s] Initial velocity of the Sun relative to the SSB at the initial epoch
    "EARTHpos_SSB_ini"  : np.array([0,0,0]), # [km]   Initial position of the Earth relative to the SSB at the initial epoch
    "EARTHvel_SSB_ini"  : np.array([0,0,0]), # [km/s] Initial velocity of the Earth relative to the SSB at the initial epoch
    "MOONpos_SSB_ini"   : np.array([0,0,0]), # [km]   Initial position of the Moon relative to the SSB at the initial epoch
    "MOONvel_SSB_ini"   : np.array([0,0,0]), # [km/s] Initial velocity of the Moon relative to the SSB at the initial epoch
    "MARSpos_SSB_ini"   : np.array([0,0,0]), # [km]   Initial position of the Mars relative to the SSB at the initial epoch
    "MARSvel_SSB_ini"   : np.array([0,0,0]), # [km/s] Initial velocity of the Mars relative to the SSB at the initial epoch
    "DEIMOSpos_SSB_ini" : np.array([0,0,0]), # [km]   Initial position of the Deimos relative to the SSB at the initial epoch
    "DEIMOSvel_SSB_ini" : np.array([0,0,0]), # [km/s] Initial velocity of the Deimos relative to the SSB at the initial epoch
    "PHOBOSpos_SSB_ini" : np.array([0,0,0]), # [km]   Initial position of the Phobos relative to the SSB at the initial epoch
    "PHOBOSvel_SSB_ini" : np.array([0,0,0]), # [km/s] Initial velocity of the Phobos relative to the SSB at the initial epoch
}
