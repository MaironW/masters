import numpy as np

# Parameters for Module DYN_MARS

DYN_MARS_par = {
    "MARSpos_SSB_ini"   : np.array([0,0,0]), # [km] Initial position of the Mars relative to the SSB at the initial epoch
    "MARSvel_SSB_ini"   : np.array([0,0,0]), # [km] Initial velocity of the Mars relative to the SSB at the initial epoch
    "DEIMOSpos_SSB_ini" : np.array([0,0,0]), # [km] Initial position of the Deimos relative to the SSB at the initial epoch
    "DEIMOSvel_SSB_ini" : np.array([0,0,0]), # [km] Initial velocity of the Deimos relative to the SSB at the initial epoch
    "PHOBOSpos_SSB_ini" : np.array([0,0,0]), # [km] Initial position of the Phobos relative to the SSB at the initial epoch
    "PHOBOSvel_SSB_ini" : np.array([0,0,0]), # [km] Initial velocity of the Phobos relative to the SSB at the initial epoch
    "MCIq_MAR_ini"      : np.array([0,0,0,0]), # Initial quaternion from MAR to MCI at the initial epoch
}