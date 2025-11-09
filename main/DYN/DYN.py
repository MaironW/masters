# Level 1 Module DYN
# Simulates the dynamic behavior of the spacecraft
# Inputs: DYN initial conditions and its own outputs (for propagation)

from .DYN_TIME  import DYN_TIME
from .DYN_SUN   import DYN_SUN
from .DYN_EARTH import DYN_EARTH
from .DYN_MARS  import DYN_MARS
from .DYN_ATT   import DYN_ATT
from .DYN_GRV   import DYN_GRV
from .DYN_TRA   import DYN_TRA

# Module output dictionary

def initialize():
    DYN_out = {
        "DYN_TIME"  : DYN_TIME.initialize(),
        "DYN_SUN"   : DYN_SUN.initialize(),
        "DYN_EARTH" : DYN_EARTH.initialize(),
        "DYN_MARS"  : DYN_MARS.initialize(),
        "DYN_ATT"   : DYN_ATT.initialize(),
        "DYN_GRV"   : DYN_GRV.initialize(),
        "DYN_TRA"   : DYN_TRA.initialize(),
    }
    return DYN_out

# Update time-dependent, non-integrated Level-2 modules
def update_algebraic(t, DYN_out, inputs):
    DYN_out = DYN_TIME.outputs(t, DYN_out)
    DYN_out = DYN_SUN.outputs(t, DYN_out)
    DYN_out = DYN_EARTH.outputs(t, DYN_out)
    DYN_out = DYN_MARS.outputs(t, DYN_out)
    DYN_out = DYN_ATT.outputs(t, DYN_out, inputs)
    DYN_out = DYN_TRA.outputs(t, DYN_out)
    DYN_out = DYN_GRV.outputs(t, DYN_out)
    return DYN_out

# Return a dict of modules that have dynamic (integrated) states
def get_dynamic_modules():
    return {
        "DYN_TIME"  : DYN_TIME,
        "DYN_SUN"   : DYN_SUN,
        "DYN_EARTH" : DYN_EARTH,
        "DYN_MARS"  : DYN_MARS,
        "DYN_TRA"   : DYN_TRA,
        "DYN_GRV"   : DYN_GRV,
    }
