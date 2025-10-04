# Level 1 Module DYN
# Simulates the dynamic behavior of the spacecraft
# Inputs: DYN initial conditions and its own outputs (for propagation)

from .DYN_TIME  import DYN_TIME
from .DYN_SUN   import DYN_SUN
from .DYN_EARTH import DYN_EARTH
from .DYN_MARS  import DYN_MARS
from .DYN_ATT   import DYN_ATT
from .DYN_TRA   import DYN_TRA

# Module output dictionary
DYN_out = {
    "DYN_TIME"  : DYN_TIME.DYN_TIME_out,
    "DYN_SUN"   : DYN_SUN.DYN_SUN_out,
    "DYN_EARTH" : DYN_EARTH.DYN_EARTH_out,
    "DYN_MARS"  : DYN_MARS.DYN_MARS_out,
    "DYN_ATT"   : DYN_ATT.DYN_ATT_out,
    "DYN_TRA"   : DYN_TRA.DYN_TRA_out,
}

# Module main function
def run():
    DYN_TIME_out  = DYN_TIME.run()
    DYN_SUN_out   = DYN_SUN.run(DYN_TIME_out)
    DYN_EARTH_out = DYN_EARTH.run(DYN_TIME_out)
    DYN_MARS_out  = DYN_MARS.run(DYN_TIME_out)
    DYN_ATT_out   = DYN_ATT.run()
    DYN_TRA_out   = DYN_TRA.run()

    # Attribute outputs to DYN output
    DYN_out["DYN_TIME"]  = DYN_TIME_out
    DYN_out["DYN_SUN"]   = DYN_SUN_out
    DYN_out["DYN_EARTH"] = DYN_EARTH_out
    DYN_out["DYN_MARS"]  = DYN_MARS_out
    DYN_out["DYN_ATT"]   = DYN_ATT_out
    DYN_out["DYN_TRA"]   = DYN_TRA_out
    return dict(DYN_out)
