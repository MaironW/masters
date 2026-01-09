# Level 1 Module DYN
# Simulates the dynamic behavior of the spacecraft
# Inputs: DYN initial conditions and its own outputs (for propagation)

from Utils.level1module   import Level1Module
from .DYN_TIME.DYN_TIME   import DYN_TIME
from .DYN_SUN.DYN_SUN     import DYN_SUN
from .DYN_EARTH.DYN_EARTH import DYN_EARTH
from .DYN_MARS.DYN_MARS   import DYN_MARS
from .DYN_ATT.DYN_ATT     import DYN_ATT
from .DYN_GRV.DYN_GRV     import DYN_GRV
from .DYN_TRA.DYN_TRA     import DYN_TRA
from .DYN_STR.DYN_STR     import DYN_STR

class DYN(Level1Module):
    def __init__(self, par_override=None):

        par_override = par_override or {}

        # Instantiate modules
        self.DYN_TIME  = DYN_TIME(par_override.get("DYN_TIME"))
        self.DYN_SUN   = DYN_SUN(par_override.get("DYN_SUN"))
        self.DYN_EARTH = DYN_EARTH(par_override.get("DYN_EARTH"))
        self.DYN_MARS  = DYN_MARS(par_override.get("DYN_MARS"))
        self.DYN_ATT   = DYN_ATT(par_override.get("DYN_ATT"))
        self.DYN_TRA   = DYN_TRA(par_override.get("DYN_TRA"))
        self.DYN_GRV   = DYN_GRV(par_override.get("DYN_GRV"))
        self.DYN_STR   = DYN_STR(par_override.get("DYN_STR"))

        # Register modules
        self.modules = [
            self.DYN_TIME,
            self.DYN_SUN,
            self.DYN_EARTH,
            self.DYN_MARS,
            self.DYN_ATT,
            self.DYN_TRA,
            self.DYN_GRV,
            self.DYN_STR,
        ]

        # Initialize all modules
        for m in self.modules:
            m.initialize(self.snapshot())

        print("DYN Module Initialized.")
