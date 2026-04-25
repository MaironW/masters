# Level 2 Module NAV_EKF
# Estimate state with Extended Kalman Filter loading NAV_CEL, NAV_PSR and NAV_CMB data

import copy
import numpy as np

from Utils.constants import CONSTANTS_par
from Utils.level2module import Level2Module
from .NAV_EKF_par import NAV_EKF_par

class NAV_EKF(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(NAV_EKF_par)
        # Apply user overrides
        if par_override is not None:
            par.update(par_override)
        # Set initial dummy state
        self.state = {
            "x_est" : None, # (n,)
            "P"     : None, # (n, n)
        }
        # Track when the last update occured
        self.last_update_time = {
            "NAV_CEL" : -np.inf,
            "NAV_PSR" : -np.inf,
            "NAV_CMB" : -np.inf,
        }
        super().__init__("NAV_EKF", par)

    # Initialization
    def initialize(self, SEN_states, NAV_states):
        # Let NAV_EKF be integrated
        self.is_dynamic = True

        self.state["x_est"] = self.par["x_est_ini"]
        self.state["P"]     = self.par["P_ini"]

        # Store the innovation per navigation method
        self.state["y_NAV_CEL"] = self.par["y_ini"]
        self.state["y_NAV_PSR"] = self.par["y_ini"]
        self.state["y_NAV_CMB"] = self.par["y_ini"]

        return self.state

    # Module main function
    # Equivalent to the update/correction step
    def update_algebraic(self, t, SEN_states, NAV_states, inputs=None):

        navigation_methods = [
            {
                "name" : "NAV_CEL",
                "data" : NAV_states["NAV_CEL"],
            },
            {
                "name" : "NAV_PSR",
                "data" : NAV_states["NAV_PSR"],
            },
            {
                "name" : "NAV_CMB",
                "data" : NAV_states["NAV_CMB"],
            },
        ]

        x_est  = self.state["x_est"]
        Pxx    = self.state["P"]

        for navigation_method in navigation_methods:
            name = navigation_method["name"]
            data = navigation_method["data"]

            # Skip if outflg is not valid
            if data[name + "outflg"] != 1:
                continue

            # Check if measurement is new
            t_valid = data["time_valid"]
            if t_valid <= self.last_update_time[name]:
                continue

            # Get measurement and measurement covariance
            z = data["z"]
            R = data["R"]

            # Check for valid measurements
            valid_measurements = ~np.isnan(z)

            # Skip if measurements are NaN (just in case)
            if np.all(~valid_measurements):
                continue

            # Get sensor-specific model
            h_fun = data["h"]
            H_fun = data["H"]

            # EKF Update
            H     = H_fun(x_est)
            z_est = h_fun(x_est)

            # Filter by measurements
            z     = z[valid_measurements]
            z_est = z_est[valid_measurements]
            H     = H[valid_measurements, :]
            R     = R[np.ix_(valid_measurements,valid_measurements)]

            # Compute inovation
            y = z - z_est

            # Compute gain
            Pzz  = H @ Pxx @ H.T + R         # Innovation covariance
            Pzz  = 0.5 * (Pzz + Pzz.T)
            Pzz += 1e-12 * np.eye(Pzz.shape[0])
            Pxy  = Pxx @ H.T                 # Cross-covariance state-measurement
            K    = np.linalg.solve(Pzz.T, Pxy.T).T # Kalman gain

            # Compute normalized innovation squared
            nis = y.T @ np.linalg.solve(Pzz, y)
            self.state[f"y_{name}"] = nis

            # Estimate
            x_est = x_est + K @ y
            I     = np.eye(Pxx.shape[0])
            Pxx   = (I - K @ H) @ Pxx @ (I - K @ H).T + K @ R @ K.T
            Pxx   = 0.5 * (Pxx + Pxx.T)
            Pxx  += 1e-12 * np.eye(Pxx.shape[0])

            # Save last update time
            self.last_update_time[name] = t_valid

        # Save final state
        self.state["x_est"] = x_est
        self.state["P"]     = Pxx
        return self.state

    # Module computation of derivatives to be integrated
    # Equivalent to the prediction step
    def derivatives(self, t, states):
        x_est = self.state["x_est"]
        Pxx   = self.state["P"]

        # Control input is none, but kept to maintain the filter structure
        u = None

        # State dynamics
        f_fun = self.f
        F_fun = self.F
        x_dot = f_fun(x_est, u, states)

        # Jacobian
        F = F_fun(x_est, states)

        # Covariance dynamics
        G = self.par["G"]
        Q = self.par["Q"]

        P_dot = F @ Pxx + Pxx @ F.T + G @ Q @ G.T

        return np.hstack([x_dot, P_dot.flatten()])

    # Return integrated variables
    def get_state(self):
        return np.hstack((self.state["x_est"], self.state["P"].flatten()))

    # Update integrated variables into the state dict
    def set_state(self, vec):
        n_states = self.par["n_states"]
        x_est = vec[:n_states]
        P     = vec[n_states:].reshape((n_states, n_states))
        P     = 0.5 * (P + P.T)
        P    += 1e-12 * np.eye(P.shape[0])
        self.state["x_est"] = x_est
        self.state["P"]     = P
        return self.state

    # Dynamics model
    def f(self, x, u, NAV_states):
        SCpos_SSB = x[0:3]
        SCvel_SSB = x[3:6]

        # Spacecraft position relative to bodies
        SCpos_SCI = SCpos_SSB - NAV_states["NAV_EPH"]["SUNpos_SSB"]   # [km]
        SCpos_ECI = SCpos_SSB - NAV_states["NAV_EPH"]["EARTHpos_SSB"] # [km]
        SCpos_MCI = SCpos_SSB - NAV_states["NAV_EPH"]["MARSpos_SSB"]  # [km]

        # Get the standard gravitational parameter around each body
        mu_SUN_cst   = CONSTANTS_par["mu_SUN_cst"]   # [km^3/s^2]
        mu_EARTH_cst = CONSTANTS_par["mu_EARTH_cst"] # [km^3/s^2]
        mu_MARS_cst  = CONSTANTS_par["mu_MARS_cst"]  # [km^3/s^2]

        # Compute the point mass acceleration (no perturbation) in each body inertial frame
        grvacc_SUN_SCI   = -mu_SUN_cst   * SCpos_SCI/np.linalg.norm(SCpos_SCI)**3 # [km/s^2]
        grvacc_EARTH_ECI = -mu_EARTH_cst * SCpos_ECI/np.linalg.norm(SCpos_ECI)**3 # [km/s^2]
        grvacc_MARS_MCI  = -mu_MARS_cst  * SCpos_MCI/np.linalg.norm(SCpos_MCI)**3 # [km/s^2]

        # Compute the gravity acceleration in the SSB frame (add all inertial models together)
        # Check Vallado c1.4 - Barycentric form of the N-body problem (eq 1-38)
        grvacc_SSB = grvacc_SUN_SCI + grvacc_EARTH_ECI + grvacc_MARS_MCI

        # Gravity is the only source of acceleration for the spacecraft
        SCacc_SSB = grvacc_SSB

        return np.hstack([SCvel_SSB, SCacc_SSB])

    # Jacobian of the dynamics model
    def F(self, x, NAV_states):
        SCpos_SSB = x[0:3]

        # Spacecraft position relative to bodies
        SCpos_SCI = SCpos_SSB - NAV_states["NAV_EPH"]["SUNpos_SSB"]   # [km]
        SCpos_ECI = SCpos_SSB - NAV_states["NAV_EPH"]["EARTHpos_SSB"] # [km]
        SCpos_MCI = SCpos_SSB - NAV_states["NAV_EPH"]["MARSpos_SSB"]  # [km]

        # Get the standard gravitational parameter around each body
        mu_SUN_cst   = CONSTANTS_par["mu_SUN_cst"]   # [km^3/s^2]
        mu_EARTH_cst = CONSTANTS_par["mu_EARTH_cst"] # [km^3/s^2]
        mu_MARS_cst  = CONSTANTS_par["mu_MARS_cst"]  # [km^3/s^2]

        # Compute the gravity gradient
        I = np.eye(3)
        SCpos_SCI_norm = np.linalg.norm(SCpos_SCI)
        SCpos_ECI_norm = np.linalg.norm(SCpos_ECI)
        SCpos_MCI_norm = np.linalg.norm(SCpos_MCI)
        dadr_SUN   = -mu_SUN_cst   * (I / SCpos_SCI_norm**3 - 3 * np.outer(SCpos_SCI, SCpos_SCI) / SCpos_SCI_norm**5)
        dadr_EARTH = -mu_EARTH_cst * (I / SCpos_ECI_norm**3 - 3 * np.outer(SCpos_ECI, SCpos_ECI) / SCpos_ECI_norm**5)
        dadr_MARS  = -mu_MARS_cst  * (I / SCpos_MCI_norm**3 - 3 * np.outer(SCpos_MCI, SCpos_MCI) / SCpos_MCI_norm**5)

        dadr = dadr_SUN + dadr_EARTH + dadr_MARS

        F_matrix = np.zeros((6, 6))
        F_matrix[0:3, 3:6] = np.eye(3)
        F_matrix[3:6, 0:3] = dadr

        return F_matrix
