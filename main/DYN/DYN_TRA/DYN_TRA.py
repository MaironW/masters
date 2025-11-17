# Level 2 Module DYN_TRA
# Simulates the propagation of the spacecraft orbital states for the simulation

from .DYN_TRA_par import DYN_TRA_par
from Utils.constants import CONSTANTS_par
from Utils import quaternions
import numpy as np

# Module output dictionary
def initialize(DYN_EARTH_out, DYN_MARS_out, DYN_SUN_out):
    # Get initial conditions based on Keplerian elements
    BODY_ini = DYN_TRA_par["BODY_ini"]
    sma_ini  = DYN_TRA_par["sma_ini"]
    ecc_ini  = DYN_TRA_par["ecc_ini"]
    incl_ini = DYN_TRA_par["incl_ini"]
    raan_ini = DYN_TRA_par["raan_ini"]
    argp_ini = DYN_TRA_par["argp_ini"]
    tano_ini = DYN_TRA_par["tano_ini"]

    if BODY_ini == "EARTH":
        # Get respective gravitational parameter
        mu = CONSTANTS_par["mu_EARTH_cst"] # [km^3/s^2]
        # Get position and velocity in the respective body centered inertial frame
        SCpos_ECI_ini, SCvel_ECI_ini = kep2rv(sma_ini, ecc_ini, incl_ini, raan_ini, argp_ini, tano_ini, mu)
        # Convert inertial states into other inertial refernces
        SCpos_SSB_ini = SCpos_ECI_ini + DYN_EARTH_out["EARTHpos_SSB"]
        SCvel_SSB_ini = SCvel_ECI_ini + DYN_EARTH_out["EARTHvel_SSB"]
        SCpos_MCI_ini = SCpos_SSB_ini - DYN_MARS_out["MARSpos_SSB"]
        SCvel_MCI_ini = SCvel_SSB_ini - DYN_MARS_out["MARSvel_SSB"]
        SCpos_SCI_ini = SCpos_SSB_ini - DYN_SUN_out["SUNpos_SSB"]
        SCvel_SCI_ini = SCvel_SSB_ini - DYN_SUN_out["SUNvel_SSB"]
    elif BODY_ini == "MARS":
        # Get respective gravitational parameter
        mu = CONSTANTS_par["mu_MARS_cst"] # [km^3/s^2]
        # Get position and velocity in the respective body centered inertial frame
        SCpos_MCI_ini, SCvel_MCI_ini = kep2rv(sma_ini, ecc_ini, incl_ini, raan_ini, argp_ini, tano_ini, mu)
        # Convert inertial states into other inertial refernces
        SCpos_SSB_ini = SCpos_MCI_ini + DYN_MARS_out["MARSpos_SSB"]
        SCvel_SSB_ini = SCvel_MCI_ini + DYN_MARS_out["MARSvel_SSB"]
        SCpos_ECI_ini = SCpos_SSB_ini - DYN_EARTH_out["EARTHpos_SSB"]
        SCvel_ECI_ini = SCvel_SSB_ini - DYN_EARTH_out["EARTHvel_SSB"]
        SCpos_SCI_ini = SCpos_SSB_ini - DYN_SUN_out["SUNpos_SSB"]
        SCvel_SCI_ini = SCvel_SSB_ini - DYN_SUN_out["SUNvel_SSB"]
    elif BODY_ini == "SUN":
        # Get respective gravitational parameter
        mu = CONSTANTS_par["mu_SUN_cst"] # [km^3/s^2]
        # Get position and velocity in the respective body centered inertial frame
        SCpos_SCI_ini, SCvel_SCI_ini = kep2rv(sma_ini, ecc_ini, incl_ini, raan_ini, argp_ini, tano_ini, mu)
        # Convert inertial states into other inertial refernces
        SCpos_SSB_ini = SCpos_SCI_ini + DYN_SUN_out["SUNpos_SSB"]
        SCvel_SSB_ini = SCvel_SCI_ini + DYN_SUN_out["SUNvel_SSB"]
        SCpos_ECI_ini = SCpos_SSB_ini - DYN_EARTH_out["EARTHpos_SSB"]
        SCvel_ECI_ini = SCvel_SSB_ini - DYN_EARTH_out["EARTHvel_SSB"]
        SCpos_MCI_ini = SCpos_SSB_ini - DYN_MARS_out["MARSpos_SSB"]
        SCvel_MCI_ini = SCvel_SSB_ini - DYN_MARS_out["MARSvel_SSB"]

    # Get frame rotation quaternions
    TERq_ECI = DYN_EARTH_out["TERq_ECI"]
    MARq_MCI = DYN_MARS_out["MARq_MCI"]
    SUNq_SSB = DYN_SUN_out["SUNq_SSB"]

    # Convert inertial references to rotational ones
    SCpos_SUN_ini = quaternions.qvecrot(SCpos_SSB_ini, SUNq_SSB)
    SCvel_SUN_ini = quaternions.qvecrot(SCvel_SSB_ini, SUNq_SSB)
    SCpos_TER_ini = quaternions.qvecrot(SCpos_ECI_ini, TERq_ECI)
    SCvel_TER_ini = quaternions.qvecrot(SCvel_ECI_ini, TERq_ECI)
    SCpos_MAR_ini = quaternions.qvecrot(SCpos_MCI_ini, MARq_MCI)
    SCvel_MAR_ini = quaternions.qvecrot(SCvel_MCI_ini, MARq_MCI)

    DYN_TRA_out = {
        "SCpos_SSB" : SCpos_SSB_ini,
        "SCvel_SSB" : SCvel_SSB_ini,
        "SCpos_SCI" : SCpos_SCI_ini,
        "SCvel_SCI" : SCvel_SCI_ini,
        "SCpos_ECI" : SCpos_ECI_ini,
        "SCvel_ECI" : SCvel_ECI_ini,
        "SCpos_MCI" : SCpos_MCI_ini,
        "SCvel_MCI" : SCvel_MCI_ini,
        "SCpos_SUN" : SCpos_SUN_ini,
        "SCvel_SUN" : SCvel_SUN_ini,
        "SCpos_TER" : SCpos_TER_ini,
        "SCvel_TER" : SCvel_TER_ini,
        "SCpos_MAR" : SCpos_MAR_ini,
        "SCvel_MAR" : SCvel_MAR_ini,
    }
    return DYN_TRA_out

# Module main function
def outputs(t, DYN_out):
    # Get parameters and states to make code more readable
    SCpos_SSB = DYN_out["DYN_TRA"]["SCpos_SSB"]  # [km]
    SCvel_SSB = DYN_out["DYN_TRA"]["SCvel_SSB"]  # [km/s]

    # Convert SSB states into other inertial refernces
    SCpos_SCI = SCpos_SSB - DYN_out["DYN_SUN"]["SUNpos_SSB"]
    SCvel_SCI = SCvel_SSB - DYN_out["DYN_SUN"]["SUNvel_SSB"]
    SCpos_ECI = SCpos_SSB - DYN_out["DYN_EARTH"]["EARTHpos_SSB"]
    SCvel_ECI = SCvel_SSB - DYN_out["DYN_EARTH"]["EARTHvel_SSB"]
    SCpos_MCI = SCpos_SSB - DYN_out["DYN_MARS"]["MARSpos_SSB"]
    SCvel_MCI = SCvel_SSB - DYN_out["DYN_MARS"]["MARSvel_SSB"]

    # Get frame rotation quaternions
    TERq_ECI = DYN_out["DYN_EARTH"]["TERq_ECI"]
    MARq_MCI = DYN_out["DYN_MARS"]["MARq_MCI"]
    SUNq_SSB = DYN_out["DYN_SUN"]["SUNq_SSB"]

    # Convert inertial references to rotational ones
    SCpos_SUN = quaternions.qvecrot(SCpos_SSB, SUNq_SSB)
    SCvel_SUN = quaternions.qvecrot(SCvel_SSB, SUNq_SSB)
    SCpos_TER = quaternions.qvecrot(SCpos_ECI, TERq_ECI)
    SCvel_TER = quaternions.qvecrot(SCvel_ECI, TERq_ECI)
    SCpos_MAR = quaternions.qvecrot(SCpos_MCI, MARq_MCI)
    SCvel_MAR = quaternions.qvecrot(SCvel_MCI, MARq_MCI)

    DYN_out["DYN_TRA"]["SCpos_SSB"] = SCpos_SSB
    DYN_out["DYN_TRA"]["SCvel_SSB"] = SCvel_SSB
    DYN_out["DYN_TRA"]["SCpos_SCI"] = SCpos_SCI
    DYN_out["DYN_TRA"]["SCvel_SCI"] = SCvel_SCI
    DYN_out["DYN_TRA"]["SCpos_ECI"] = SCpos_ECI
    DYN_out["DYN_TRA"]["SCvel_ECI"] = SCvel_ECI
    DYN_out["DYN_TRA"]["SCpos_MCI"] = SCpos_MCI
    DYN_out["DYN_TRA"]["SCvel_MCI"] = SCvel_MCI
    DYN_out["DYN_TRA"]["SCpos_SUN"] = SCpos_SUN
    DYN_out["DYN_TRA"]["SCvel_SUN"] = SCvel_SUN
    DYN_out["DYN_TRA"]["SCpos_TER"] = SCpos_TER
    DYN_out["DYN_TRA"]["SCvel_TER"] = SCvel_TER
    DYN_out["DYN_TRA"]["SCpos_MAR"] = SCpos_MAR
    DYN_out["DYN_TRA"]["SCvel_MAR"] = SCvel_MAR

    return dict(DYN_out)

# Module computation of derivatives to be integrated
def derivatives(t, DYN_out):
    # Get parameters and states to make code more readable
    SCvel_SSB  = DYN_out["DYN_TRA"]["SCvel_SSB"]  # [km/s]
    grvacc_SSB = DYN_out["DYN_GRV"]["grvacc_SSB"] # [km/s^2]
    # Return derivatives
    dSCpos_SSB = SCvel_SSB
    dSCvel_SSB = grvacc_SSB
    return np.hstack([dSCpos_SSB, dSCvel_SSB])

# Return integrated variables
def get_state(DYN_TRA_out):
    return np.hstack((DYN_TRA_out["SCpos_SSB"], DYN_TRA_out["SCvel_SSB"]))

# Update integrated variables into the state dict
def set_state(DYN_TRA_out, vec):
    DYN_TRA_out["SCpos_SSB"] = vec[0:3]
    DYN_TRA_out["SCvel_SSB"] = vec[3:6]
    return DYN_TRA_out

# Convert Keplerian elements to cartesian position and velocity
def kep2rv(sma, ecc, incl, raan, argp, tano, mu):
    # Compute position and velocity on the perifocal frame
    p    = sma*(1 - ecc**2) # Semi-latus rectum
    r_pf = (p/(1 + ecc*np.cos(tano)))*np.array([np.cos(tano), np.sin(tano), 0.0])
    v_pf = np.sqrt(mu/p)*np.array([-np.sin(tano), ecc + np.cos(tano), 0.0])

    # Rotation from perifocal to inertial frame
    c_raan = np.cos(raan); s_raan = np.sin(raan)
    c_incl = np.cos(incl); s_incl = np.sin(incl)
    c_argp = np.cos(argp); s_argp = np.sin(argp)

    R = np.array([
        [ c_raan*c_argp - s_raan*s_argp*c_incl, -c_raan*s_argp - s_raan*c_argp*c_incl,  s_raan*s_incl],
        [ s_raan*c_argp + c_raan*s_argp*c_incl, -s_raan*s_argp + c_raan*c_argp*c_incl, -c_raan*s_incl],
        [                        s_argp*s_incl,                         c_argp*s_incl,         c_incl]
    ])

    r_ine = R.dot(r_pf)
    v_ine = R.dot(v_pf)
    return r_ine, v_ine
