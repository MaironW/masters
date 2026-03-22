# Load Pulsar data, useful to keep both DYN_PSR and SEN_PSR in sync

import numpy as np
import pandas as pd

class PulsarDatabase:
    def __init__(self, filename):
        df = pd.read_csv(filename, skipinitialspace=True)

        self.name = df["name"].values # Pulsar name
        self.lat  = df["lat"].values  # [deg] Galactic latitude
        self.lon  = df["lon"].values  # [deg] Galactic longidude
        self.f    = df["f"].values    # [Hz] Pulse frequency
        self.Fx   = df["Fx"].values   # [ph/cm^2/s] Pulsar raduation flux
        self.pf   = df["pf"].values   # Flux pulsed fraction
        self.W    = df["W"].values    # [s] Pulse width
        self.D0   = df["D0"].values   # [kpc] Pulsar distance from SSB

        self.Bx = 0.005 # [ph/cm^2/s] X-ray background radiation flux
        self.n_pulsars = len(self.name)

        self.PULSARdir_SSB = self.latlon2dir(self.lat, self.lon)

    # Convert Galatic Latitude and Longitude in degrees to a direction vector in the SSB frame
    def latlon2dir(self, lat, lon):
        lat = np.deg2rad(lat)
        lon = np.deg2rad(lon)

        clat, slat = np.cos(lat), np.sin(lat)
        clon, slon = np.cos(lon), np.sin(lon)

        dir_galatic = np.array([clat*clon, clat*slon, slat])

        # Rotation matrix from Galactic to Equatorial plane J2000
        R = np.array([
            [-0.0548755604, -0.8734370902, -0.4838350155],
            [ 0.4941094279, -0.4448296300,  0.7469822445],
            [-0.8676661490, -0.1980763734,  0.4559837762]
        ])

        return (R @ dir_galatic).T
