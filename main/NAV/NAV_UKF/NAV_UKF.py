# Level 2 Module NAV_UKF
# Estimate state with Unscented Kalman Filter loading NAV_CEL, NAV_PSR and NAV_CMB data

import copy
import numpy as np
from scipy.linalg import block_diag

from Utils.constants import CONSTANTS_par
from Utils.level2module import Level2Module
from .NAV_UKF_par import NAV_UKF_par

class NAV_UKF(Level2Module):
    def __init__(self, par_override=None):
        # Start with default parameters
        par = copy.deepcopy(NAV_UKF_par)
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
        super().__init__("NAV_UKF", par)

    # Initialization
    def initialize(self, SEN_states, NAV_states):
        # NAV_UKF should not be integrated externally
        self.is_dynamic = False

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

        # Load states and parameters
        x_est  = self.state["x_est"]
        Pxx    = self.state["P"]
        dt     = self.par["dt"]
        Q      = self.par["Q"]
        G      = self.par["G"]

        # Generate sigma points
        X_sigma_aug, Wm, Wc = self.generate_sigma_points(x_est, Pxx, Q)

        # Get individual sigma-points
        n_states  = len(x_est)
        n_process = Q.shape[0]
        X_sigma = X_sigma_aug[0 : n_states, :]
        W_sigma = X_sigma_aug[n_states : n_states + n_process, :]

        # Propagate the sigma points through integration
        X_sigma_prop = []
        for i in range(X_sigma.shape[1]):
            x_i = X_sigma[:, i]
            w_i = W_sigma[:, i]
            x_next = self.propagate_sigma(x_i, dt, NAV_states) + G @ w_i
            X_sigma_prop.append(x_next)
        X_sigma_prop = np.array(X_sigma_prop).T

        # Mean, prediction
        x_pred = X_sigma_prop @ Wm
        x_est  = np.copy(x_pred)

        # State covariance
        Pxx = np.zeros((n_states, n_states))
        for i in range(X_sigma_prop.shape[1]):
            dx = X_sigma_prop[:, i] - x_pred
            Pxx += Wc[i] * np.outer(dx, dx)
        Pxx = 0.5 * (Pxx + Pxx.T) # Symmetry fix
        Pxx += 1e-12*np.eye(n_states) # Add low value to ensure convergence

        # Sequential update iterating on the navigation methods
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

            # UKF Update
            Z_sigma_prop = []
            for i in range(X_sigma_prop.shape[1]):
                z_i = h_fun(X_sigma_prop[:, i])
                Z_sigma_prop.append(z_i)
            Z_sigma_prop = np.array(Z_sigma_prop).T
            z_est = Z_sigma_prop @ Wm

            # Filter by measurements
            z            = z[valid_measurements]
            z_est        = z_est[valid_measurements]
            Z_sigma_prop = Z_sigma_prop[valid_measurements, :]
            R      = R[np.ix_(valid_measurements,valid_measurements)]
            R      = 0.5 * (R + R.T) # Symmetry fix
            R     += 1e-12 * np.eye(R.shape[0]) # Add low value to ensure convergence
            n_mes  = len(z)

            # Covariances
            Pzz = np.zeros((n_mes,    n_mes)) # Innovation covariance
            Pxz = np.zeros((n_states, n_mes)) # Cross-covariance state-measurement

            for i in range(Z_sigma_prop.shape[1]):
                dx = X_sigma_prop[:, i] - x_est
                dz = Z_sigma_prop[:, i] - z_est

                Pzz += Wc[i] * np.outer(dz, dz)
                Pxz += Wc[i] * np.outer(dx, dz)

            Pzz += R
            Pzz += 1e-12*np.eye(n_mes) # Add low value to ensure convergence

            # Compute gain
            K = np.linalg.solve(Pzz.T, Pxz.T).T

            # Compute inovation
            y = z - z_est

            # Compute normalized innovation squared
            nis = y.T @ np.linalg.solve(Pzz, y)
            self.state[f"y_{name}"] = nis

            # Estimate
            x_est = x_est + K @ y
            Pxx   = Pxx - K @ Pzz @ K.T
            Pxx   = 0.5 *(Pxx + Pxx.T) # Symmetry fix
            Pxx  += 1e-12*np.eye(n_states) # Add low value to ensure convergence

            # Save last update time
            self.last_update_time[name] = t_valid

            # Regenerate sigma-points from last update
            X_sigma_aug, Wm, Wc = self.generate_sigma_points(x_est, Pxx, Q)
            X_sigma = X_sigma_aug[0 : n_states, :]
            W_sigma = X_sigma_aug[n_states : n_states + n_process, :]

            # Re-propagate sigma points
            X_sigma_prop = np.copy(X_sigma)

        # Save final state
        self.state["x_est"] = x_est
        self.state["P"]     = Pxx
        return self.state

    # Return integrated variables
    def get_state(self):
        return np.hstack((self.state["x_est"], self.state["P"].flatten()))

    # Update integrated variables into the state dict
    def set_state(self, vec):
        n_states = self.par["n_states"]
        x_est = vec[:n_states]
        P     = vec[n_states:].reshape((n_states, n_states))
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

    # Generate sigma-points for the UKF
    def generate_sigma_points(self, x, P, Q):
        n_states  = len(x)                # State size
        n_process = Q.shape[0]            # Process noise size
        n_a       = n_states + n_process  # Augmented state size

        # Augmented state
        x_aug = np.hstack([x, np.zeros(n_process)])

        # Augmented covariance
        P_aug = block_diag(P, Q)

        # UKF parameters following Wan & van der Merwe (2000)
        alpha  = self.par["alpha"]
        alpha2 = alpha**2
        beta   = self.par["beta"]
        kappa  = self.par["kappa"]
        L      = alpha2*(n_a + kappa) - n_a
        gamma  = np.sqrt(n_a + L)

        # Generate augmented sigma-points
        aux_sqrt = gamma * np.linalg.cholesky(P_aug + 1e-12*np.eye(n_a))

        sigma_points = [x_aug]
        for i in range(n_a):
            sigma_points.append(x_aug + aux_sqrt[:, i])
            sigma_points.append(x_aug - aux_sqrt[:, i])

        sigma_points = np.array(sigma_points).T

        # Mean weights
        Wm    = np.full(2*n_a + 1, 0.5 / (n_a + L))
        Wm[0] = L / (n_a + L)

        # Covariance weights
        Wc     = np.copy(Wm)
        Wc[0] += 1 - alpha2 + beta

        return sigma_points, Wm, Wc

    # Internal RK4 integrator for the UKF
    def propagate_sigma(self, x, dt, NAV_states):
        def f_local(x_local):
            return self.f(x_local, None, NAV_states)
        k1 = f_local(x)
        k2 = f_local(x + 0.5 * dt * k1)
        k3 = f_local(x + 0.5 * dt * k2)
        k4 = f_local(x + dt * k3)
        return x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
