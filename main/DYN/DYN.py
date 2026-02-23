# Level 1 Module DYN
# Simulates the dynamic behavior of the spacecraft
# Inputs: DYN initial conditions and its own outputs (for propagation)

from Utils.level1module import Level1Module
from .DYN_TIME.DYN_TIME import DYN_TIME
from .DYN_EPH.DYN_EPH   import DYN_EPH
from .DYN_ATT.DYN_ATT   import DYN_ATT
from .DYN_GRV.DYN_GRV   import DYN_GRV
from .DYN_TRA.DYN_TRA   import DYN_TRA
from .DYN_STR.DYN_STR   import DYN_STR
from .DYN_PSR.DYN_PSR   import DYN_PSR

class DYN(Level1Module):
    def __init__(self, par_override=None):

        par_override = par_override or {}

        # Instantiate modules
        self.DYN_TIME = DYN_TIME(par_override.get("DYN_TIME"))
        self.DYN_EPH  = DYN_EPH(par_override.get("DYN_EPH"))
        self.DYN_ATT  = DYN_ATT(par_override.get("DYN_ATT"))
        self.DYN_TRA  = DYN_TRA(par_override.get("DYN_TRA"))
        self.DYN_GRV  = DYN_GRV(par_override.get("DYN_GRV"))
        self.DYN_STR  = DYN_STR(par_override.get("DYN_STR"))
        self.DYN_PSR  = DYN_PSR(par_override.get("DYN_PSR"))

        # Register modules
        self.modules = [
            self.DYN_TIME,
            self.DYN_EPH,
            self.DYN_ATT,
            self.DYN_TRA,
            self.DYN_GRV,
            self.DYN_STR,
            self.DYN_PSR,
        ]

        # Initialize all modules
        for m in self.modules:
            m.initialize(None, self.snapshot())

        print("DYN Module Initialized.")

    # Update time-dependent, non-integrated Level-2 modules
    def update_algebraic(self, t, parent_states, inputs):
        for m in self.modules:
            m.update_algebraic(t, None, self.snapshot(), inputs)
