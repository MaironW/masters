# Load Pulsar data, useful to keep both DYN_PSR and SEN_PSR in sync

import numpy as np
import pandas as pd

class PulsarDatabase:
    def __init__(self, filename):
        df = pd.read_csv(filename)

        self.name  = df["name"].values #  Pulsar name
        self.ra    = df["ra"].values    # [deg] Right ascension from SSB
        self.dec   = df["dec"].values   # [deg] Declination from SSB
        self.epoch = df["epoch"].values # [MJD] Epoch for frequency
        self.f     = df["f"].values     # [Hz] Pulse frequency
        self.df    = df["df"].values    # [Hz/s] First derivative of pulse
        self.Fx    = df["Fx"].values    # [photons/m^2/s] Pulsar raduation flux
        self.pf    = df["pf"].values    # Flux pulsed fraction
        self.d     = df["d"].values     # Pulse duty cycle
        self.D0    = df["D0"].values    # [kpc] Pulsar distance from SSB

        self.Bx = 50 # [photon/m^2/s] X-ray background radiation flux
        self.n_pulsars = len(self.name)

        self.PULSARdir_SSB = self.radec2dir(self.ra, self.dec)

    # Convert Right Ascension and Declination in degrees to a direction vector
    def radec2dir(self, ra, dec):
        ra  = np.deg2rad(ra)
        dec = np.deg2rad(dec)
        x = np.cos(dec)*np.cos(ra)
        y = np.cos(dec)*np.sin(ra)
        z = np.sin(dec)
        return np.vstack([x,y,z]).T
