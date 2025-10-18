import numpy as np

# Parameters for Module DYN_GRV

DYN_GRV_par = {
    "grvacc_TER_ini" : np.array([0,0,0]),  # [km/s^2] Initial gravity acceleration of the spacecraft around Earth on TER frame
    "grvacc_ECI_ini" : np.array([0,0,0]),  # [km/s^2] Initial gravity acceleration of the spacecraft around Earth on ECI frame
    "grvacc_MAR_ini" : np.array([0,0,0]),  # [km/s^2] Initial gravity acceleration of the spacecraft around Mars on MAR frame
    "grvacc_MCI_ini" : np.array([0,0,0]),  # [km/s^2] Initial gravity acceleration of the spacecraft around Mars on MCI frame
    "grvacc_SUN_ini" : np.array([0,0,0]),  # [km/s^2] Initial gravity acceleration of the spacecraft around the Sun on SUN frame
    "grvacc_SSB_ini" : np.array([0,0,0])   # [km/s^2] Initial gravity acceleration of the spacecraft around the Sun on SSB frame
}
