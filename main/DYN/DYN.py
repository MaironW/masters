# Level 1 Module DYN
# Simulates the dynamic behavior of the spacecraft
# Inputs: DYN initial conditions and its own outputs (for propagation)

from .DYN_TIME.DYN_TIME   import DYN_TIME
from .DYN_SUN.DYN_SUN     import DYN_SUN
from .DYN_EARTH.DYN_EARTH import DYN_EARTH
# from .DYN_MARS.DYN_MARS   import DYN_MARS
# from .DYN_ATT.DYN_ATT     import DYN_ATT
# from .DYN_GRV.DYN_GRV     import DYN_GRV
# from .DYN_TRA.DYN_TRA     import DYN_TRA

class DYN:
    def __init__(self, par_override=None):

        par_override = par_override or {}

        # Instantiate modules
        self.DYN_TIME  = DYN_TIME(par_override.get("DYN_TIME"))
        self.DYN_SUN   = DYN_SUN(par_override.get("DYN_SUN"))
        self.DYN_EARTH = DYN_EARTH(par_override.get("DYN_EARTH"))

        # Register modules
        self.modules = [
            self.DYN_TIME,
            self.DYN_SUN,
            self.DYN_EARTH
        ]
   
        # Initialize all modules
        for m in self.modules:
            m.initialize(self.snapshot())

        print("DYN Module Initialized.")

    # Organize dynamic modules
    @property
    def dynamic_modules(self):
        return [m for m in self.modules if m.is_dynamic]

    # Update time-dependent, non-integrated Level-2 modules
    def update_algebraic(self, t):
        for m in self.modules:
            m.update_algebraic(t, self.snapshot())

    # Return a dict of modules that have dynamic (integrated) states
    def get_dynamic_modules(self):
        return self.dynamic_modules

    def get_state(self):
        return [m.get_state() for m in self.dynamic_modules]

    def set_state(self, vecs):
        for m, v in zip(self.dynamic_modules, vecs):
            m.set_state(v)

    def derivatives(self, t):
        return [m.derivatives(t) for m in self.dynamic_modules]

    def snapshot(self):
        return {m.name: m.state.copy() for m in self.modules}