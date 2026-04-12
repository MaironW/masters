# Load Pulsar data, useful to keep both DYN_PSR and SEN_PSR in sync

import pandas as pd
from Utils import misc

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

        self.PULSARdir_SSB = misc.GALtoSSB(misc.latlon2dir(self.lat, self.lon))
