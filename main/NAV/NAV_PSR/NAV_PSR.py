# Level 2 Module NAV_PSR
# Provides the inputs for a Kalman Filter based on Pulsar measurements
# Inputs: Selected pulsars TOAs difference from the SC to the SSB

import copy
import numpy as np

from Utils.constants import CONSTANTS_par
from Utils.level2module import Level2Module
from Utils.pulsar_database import PulsarDatabase
from .NAV_PSR_par import NAV_PSR_par

class NAV_PSR(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(NAV_PSR_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "NAV_PSRoutflg" : par["NAV_PSRoutflg_ini"],
        }
        super().__init__("NAV_PSR", par)

    # Initialization
    def initialize(self, SEN_states, NAV_states):
         # Load Pulsar database
        kernel_dir = "Utils/kernels/"
        pulsar_data_file = kernel_dir + "pulsar.csv"
        pulsar_data = PulsarDatabase(pulsar_data_file)

        # Pulsar parameters should not change over the simulation
        self.par["name"]           = pulsar_data.name          # Pulsar name
        self.par["n_pulsars"]      = pulsar_data.n_pulsars     # Number of pulsars
        self.par["D0"]             = pulsar_data.D0*CONSTANTS_par["pc2m_cst"] # [km] Pulsar distance from SSB
        self.par["PULSARSdir_SSB"] = pulsar_data.PULSARdir_SSB # Direction of Pulsar from SSB

        # Allocate initial values
        n_pulsars = self.par["n_pulsars"]
        self.par["PULSARSid_mes_ini"] = np.full(n_pulsars, None)
        self.par["z_ini"]             = np.full(n_pulsars, np.nan)
        self.par["R_ini"]             = np.full((n_pulsars, n_pulsars), np.nan)

        self.state["PULSARSid_mes"]  = self.par["PULSARSid_mes_ini"]
        self.state["z"]              = self.par["z_ini"]
        self.state["R"]              = self.par["R_ini"]
        self.state["SSBpos_SUN_ref"] = self.par["SSBpos_SUN_ref_ini"]

        # Update initial state
        self.state = self.update_algebraic(0, SEN_states, NAV_states)

        return self.state

    # Module main function
    def update_algebraic(self, t, SEN_states, NAV_states, inputs=None):
        SEN_PSRoutflg  = SEN_states["SEN_PSR"]["SEN_PSRoutflg"]
        SUNpos_SSB_ref = NAV_states["NAV_EPH"]["SUNpos_SSB"]

        # Only update outputs if PSRoutflg is valid
        # Otherwise, return default values
        if SEN_PSRoutflg == 0:
            self.state["NAV_PSRoutflg"]  = self.par["NAV_PSRoutflg_ini"]
            self.state["z"]              = self.par["z_ini"]
            self.state["R"]              = self.par["R_ini"]
            self.state["SSBpos_SUN_ref"] = self.par["SSBpos_SUN_ref_ini"]

        # PSR output is valid
        else:
            # Compute the position of SSB relative to the Sun
            SSBpos_SUN_ref = -SUNpos_SSB_ref

            # Load sensor outputs
            OBTdt_TDB_mes  = SEN_states["SEN_PSR"]["OBTdt_TDB_mes"]

            n_pulsars = self.par["n_pulsars"]
            z = np.full(n_pulsars, np.nan)
            R = np.full(n_pulsars, np.nan)

            # Compute dt if pulsars are visible
            sigma_TOA = np.asarray(self.par["sigma_TOA"])
            z = OBTdt_TDB_mes
            R = sigma_TOA**2

            if np.sum(~np.isnan(z)) >= 3:
                z = z
                R = np.diag(R)
                NAV_PSRoutflg = 1
            else:
                z = self.par["z_ini"]
                R = np.diag(self.par["R_ini"])
                NAV_PSRoutflg = 0

            # Update states
            self.state["NAV_PSRoutflg"]  = NAV_PSRoutflg
            self.state["z"]              = z
            self.state["R"]              = R
            self.state["SSBpos_SUN_ref"] = SSBpos_SUN_ref

        return self.state

    # Return the pulsar x-ray navigation measurement model for Kalman filtering
    def h(self, x):
        # x = [SCpos_SSB, SCvel_SSB]
        SCpos_SSB = x[0:3]

        # Load parameters and states
        mu_SUN_cst      = CONSTANTS_par["mu_SUN_cst"] # [km^3/s^2]
        light_speed_cst = CONSTANTS_par["light_speed_cst"] # [km/s]
        SSBpos_SUN      = self.state["SSBpos_SUN_ref"] # [km]
        PULSARSdir_SSB  = self.par["PULSARSdir_SSB"]
        D0              = self.par["D0"] # [km]

        # Compute the time of arrival as perceived by the spacecraft
        n_dot_r = PULSARSdir_SSB @ SCpos_SSB
        r_dot_r = SCpos_SSB @ SCpos_SSB
        n_dot_b = PULSARSdir_SSB @ SSBpos_SUN
        b_dot_r = SSBpos_SUN @ SCpos_SSB

        # Roemer delay
        doppler_delay = n_dot_r / light_speed_cst # [s]
        annual_parallax_delay = 1/(2*light_speed_cst*D0) * (n_dot_r**2 - r_dot_r + 2*n_dot_b*n_dot_r - 2*b_dot_r) # [s]
        roemer_delay = doppler_delay + annual_parallax_delay # [s]

        # Shapiro delay considering only the effect of the Sun as the major source of spacetime curvature in the Solar System
        SCpos_SSB_norm  = np.linalg.norm(SCpos_SSB)
        SSBpos_SUN_norm = np.linalg.norm(SSBpos_SUN)
        num = n_dot_r + SCpos_SSB_norm
        den = n_dot_b + SSBpos_SUN_norm
        shapiro_delay = 2*mu_SUN_cst/light_speed_cst**3 * np.log(np.abs(num/den + 1)) # [s]

        h_vec = roemer_delay + shapiro_delay

        return h_vec

    # Return the Jacobian of the measurment model for EKF
    def H(self, x):
        # x = [SCpos_SSB, SCvel_SSB]
        SCpos_SSB = x[0:3]

        # Load parameters and states
        mu_SUN_cst      = CONSTANTS_par["mu_SUN_cst"] # [km^3/s^2]
        light_speed_cst = CONSTANTS_par["light_speed_cst"] # [km/s]
        SSBpos_SUN      = self.state["SSBpos_SUN_ref"] # [km]
        PULSARSdir_SSB  = self.par["PULSARSdir_SSB"]
        D0              = self.par["D0"] # [km]
        n_pulsars       = self.par["n_pulsars"]
        n_states        = len(x)

        H_matrix = np.zeros((n_pulsars, n_states))

        # Compute scalars
        SCpos_SSB_norm  = np.linalg.norm(SCpos_SSB)
        SSBpos_SUN_norm = np.linalg.norm(SSBpos_SUN)

        # Compute the time of arrival as perceived by the spacecraft
        n_dot_r = PULSARSdir_SSB @ SCpos_SSB
        n_dot_b = PULSARSdir_SSB @ SSBpos_SUN

        for i in range(n_pulsars):

            # Doppler delay
            dh1dr = PULSARSdir_SSB[i]/light_speed_cst

            # Annual parallax delay
            dh2dr = 1/(light_speed_cst*D0[i]) * (n_dot_r[i]*PULSARSdir_SSB[i] - SCpos_SSB + n_dot_b[i]*PULSARSdir_SSB[i] - SSBpos_SUN)

            # Shappiro delay
            dh3dr = 2*mu_SUN_cst/light_speed_cst**3 / (n_dot_r[i] + SCpos_SSB_norm + n_dot_b[i] + SSBpos_SUN_norm) * (PULSARSdir_SSB[i] + SCpos_SSB/SCpos_SSB_norm)

            # Fill matrix
            H_matrix[i, 0:3] = dh1dr + dh2dr + dh3dr

        return H_matrix
