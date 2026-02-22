# Level 2 Module DYN_TRA
# Simulates the propagation of the spacecraft orbital states for the simulation

import copy
import numpy as np

from Utils.constants import CONSTANTS_par
from Utils.level2module import Level2Module
from .DYN_TRA_par import DYN_TRA_par

class DYN_TRA(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(DYN_TRA_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "SCpos_SSB" : par["SCpos_SSB_ini"],
            "SCvel_SSB" : par["SCvel_SSB_ini"],
        }
        super().__init__("DYN_TRA", par)

    # Initialization
    def initialize(self, DYN_states):
        # Let DYN_TRA be integrated
        self.is_dynamic = True

        if self.par["ref_elements"] == "kep":
            self.state = self.initialize_keplerian(DYN_states)
        else: # "rvi"
            self.state = self.initialize_rvi(DYN_states)

        return self.state

    # Module main function
    def update_algebraic(self, t, DYN_states, inputs=None):
        # Get parameters and states to make code more readable
        SCpos_SSB = DYN_states["DYN_TRA"]["SCpos_SSB"]  # [km]
        SCvel_SSB = DYN_states["DYN_TRA"]["SCvel_SSB"]  # [km/s]

        # Convert SSB states into other inertial refernces
        SCpos_SCI = SCpos_SSB - DYN_states["DYN_EPH"]["SUNpos_SSB"]
        SCvel_SCI = SCvel_SSB - DYN_states["DYN_EPH"]["SUNvel_SSB"]
        SCpos_ECI = SCpos_SSB - DYN_states["DYN_EPH"]["EARTHpos_SSB"]
        SCvel_ECI = SCvel_SSB - DYN_states["DYN_EPH"]["EARTHvel_SSB"]
        SCpos_MCI = SCpos_SSB - DYN_states["DYN_EPH"]["MARSpos_SSB"]
        SCvel_MCI = SCvel_SSB - DYN_states["DYN_EPH"]["MARSvel_SSB"]

        self.state["SCpos_SSB"] = SCpos_SSB
        self.state["SCvel_SSB"] = SCvel_SSB
        self.state["SCpos_SCI"] = SCpos_SCI
        self.state["SCvel_SCI"] = SCvel_SCI
        self.state["SCpos_ECI"] = SCpos_ECI
        self.state["SCvel_ECI"] = SCvel_ECI
        self.state["SCpos_MCI"] = SCpos_MCI
        self.state["SCvel_MCI"] = SCvel_MCI

        return self.state

    # Module computation of derivatives to be integrated
    def derivatives(self, t, DYN_states):
        # Get parameters and states to make code more readable
        SCvel_SSB  = self.state["SCvel_SSB"]             # [km/s]
        grvacc_SSB = DYN_states["DYN_GRV"]["grvacc_SSB"] # [km/s^2]
        # Return derivatives
        dSCpos_SSB = SCvel_SSB
        dSCvel_SSB = grvacc_SSB
        return np.hstack([dSCpos_SSB, dSCvel_SSB])

    # Return integrated variables
    def get_state(self):
        return np.hstack((self.state["SCpos_SSB"], self.state["SCvel_SSB"]))

    # Update integrated variables into the state dict
    def set_state(self, vec):
        self.state["SCpos_SSB"] = vec[0:3]
        self.state["SCvel_SSB"] = vec[3:6]
        return self.state

    # Use Keplerian Elements to initialize DYN_TRA
    def initialize_keplerian(self, DYN_states):
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
            SCpos_ECI_ini, SCvel_ECI_ini = self.kep2rvi(sma_ini, ecc_ini, incl_ini, raan_ini, argp_ini, tano_ini, mu)
            # Convert inertial states into other inertial refernces
            SCpos_SSB_ini = SCpos_ECI_ini + DYN_states["DYN_EPH"]["EARTHpos_SSB"]
            SCvel_SSB_ini = SCvel_ECI_ini + DYN_states["DYN_EPH"]["EARTHvel_SSB"]
            SCpos_MCI_ini = SCpos_SSB_ini - DYN_states["DYN_EPH"]["MARSpos_SSB"]
            SCvel_MCI_ini = SCvel_SSB_ini - DYN_states["DYN_EPH"]["MARSvel_SSB"]
            SCpos_SCI_ini = SCpos_SSB_ini - DYN_states["DYN_EPH"]["SUNpos_SSB"]
            SCvel_SCI_ini = SCvel_SSB_ini - DYN_states["DYN_EPH"]["SUNvel_SSB"]
        elif BODY_ini == "MARS":
            # Get respective gravitational parameter
            mu = CONSTANTS_par["mu_MARS_cst"] # [km^3/s^2]
            # Get position and velocity in the respective body centered inertial frame
            SCpos_MCI_ini, SCvel_MCI_ini = self.kep2rvi(sma_ini, ecc_ini, incl_ini, raan_ini, argp_ini, tano_ini, mu)
            # Convert inertial states into other inertial refernces
            SCpos_SSB_ini = SCpos_MCI_ini + DYN_states["DYN_EPH"]["MARSpos_SSB"]
            SCvel_SSB_ini = SCvel_MCI_ini + DYN_states["DYN_EPH"]["MARSvel_SSB"]
            SCpos_ECI_ini = SCpos_SSB_ini - DYN_states["DYN_EPH"]["EARTHpos_SSB"]
            SCvel_ECI_ini = SCvel_SSB_ini - DYN_states["DYN_EPH"]["EARTHvel_SSB"]
            SCpos_SCI_ini = SCpos_SSB_ini - DYN_states["DYN_EPH"]["SUNpos_SSB"]
            SCvel_SCI_ini = SCvel_SSB_ini - DYN_states["DYN_EPH"]["SUNvel_SSB"]
        elif BODY_ini == "SUN":
            # Get respective gravitational parameter
            mu = CONSTANTS_par["mu_SUN_cst"] # [km^3/s^2]
            # Get position and velocity in the respective body centered inertial frame
            SCpos_SCI_ini, SCvel_SCI_ini = self.kep2rvi(sma_ini, ecc_ini, incl_ini, raan_ini, argp_ini, tano_ini, mu)
            # Convert inertial states into other inertial refernces
            SCpos_SSB_ini = SCpos_SCI_ini + DYN_states["DYN_EPH"]["SUNpos_SSB"]
            SCvel_SSB_ini = SCvel_SCI_ini + DYN_states["DYN_EPH"]["SUNvel_SSB"]
            SCpos_ECI_ini = SCpos_SSB_ini - DYN_states["DYN_EPH"]["EARTHpos_SSB"]
            SCvel_ECI_ini = SCvel_SSB_ini - DYN_states["DYN_EPH"]["EARTHvel_SSB"]
            SCpos_MCI_ini = SCpos_SSB_ini - DYN_states["DYN_EPH"]["MARSpos_SSB"]
            SCvel_MCI_ini = SCvel_SSB_ini - DYN_states["DYN_EPH"]["MARSvel_SSB"]

        # Update output
        self.state = {
            "SCpos_SSB" : SCpos_SSB_ini,
            "SCvel_SSB" : SCvel_SSB_ini,
            "SCpos_ECI" : SCpos_ECI_ini,
            "SCvel_ECI" : SCvel_ECI_ini,
            "SCpos_MCI" : SCpos_MCI_ini,
            "SCvel_MCI" : SCvel_MCI_ini,
            "SCpos_SCI" : SCpos_SCI_ini,
            "SCvel_SCI" : SCvel_SCI_ini,
        }

        return self.state

    # Use state vectors (position and velocitu) in the inertial frame to initialize DYN_TRA
    def initialize_rvi(self, DYN_states):
        SCpos_SSB_ini = self.par["SCpos_SSB_ini"]
        SCvel_SSB_ini = self.par["SCvel_SSB_ini"]
        SCpos_ECI_ini = SCpos_SSB_ini - DYN_states["DYN_EPH"]["EARTHpos_SSB"]
        SCvel_ECI_ini = SCvel_SSB_ini - DYN_states["DYN_EPH"]["EARTHvel_SSB"]
        SCpos_MCI_ini = SCpos_SSB_ini - DYN_states["DYN_EPH"]["MARSpos_SSB"]
        SCvel_MCI_ini = SCvel_SSB_ini - DYN_states["DYN_EPH"]["MARSvel_SSB"]
        SCpos_SCI_ini = SCpos_SSB_ini - DYN_states["DYN_EPH"]["SUNpos_SSB"]
        SCvel_SCI_ini = SCvel_SSB_ini - DYN_states["DYN_EPH"]["SUNvel_SSB"]

        # Update output
        self.state = {
            "SCpos_SSB" : SCpos_SSB_ini,
            "SCvel_SSB" : SCvel_SSB_ini,
            "SCpos_ECI" : SCpos_ECI_ini,
            "SCvel_ECI" : SCvel_ECI_ini,
            "SCpos_MCI" : SCpos_MCI_ini,
            "SCvel_MCI" : SCvel_MCI_ini,
            "SCpos_SCI" : SCpos_SCI_ini,
            "SCvel_SCI" : SCvel_SCI_ini,
        }

        return self.state

    # Convert Keplerian elements to cartesian position and velocity
    def kep2rvi(self, sma, ecc, incl, raan, argp, tano, mu):
        # Compute position and velocity on the perifocal frame
        p     = sma*(1 - ecc**2) # Semi-latus rectum
        r_PQW = (p/(1 + ecc*np.cos(tano)))*np.array([np.cos(tano), np.sin(tano), 0.0])
        v_PQW = np.sqrt(mu/p)*np.array([-np.sin(tano), ecc + np.cos(tano), 0.0])

        # Rotation from perifocal to inertial frame
        c_raan = np.cos(raan); s_raan = np.sin(raan)
        c_incl = np.cos(incl); s_incl = np.sin(incl)
        c_argp = np.cos(argp); s_argp = np.sin(argp)

        R = np.array([
            [ c_raan*c_argp - s_raan*s_argp*c_incl, -c_raan*s_argp - s_raan*c_argp*c_incl,  s_raan*s_incl],
            [ s_raan*c_argp + c_raan*s_argp*c_incl, -s_raan*s_argp + c_raan*c_argp*c_incl, -c_raan*s_incl],
            [                        s_argp*s_incl,                         c_argp*s_incl,         c_incl]
        ])

        r_ine = R.dot(r_PQW)
        v_ine = R.dot(v_PQW)
        return r_ine, v_ine
