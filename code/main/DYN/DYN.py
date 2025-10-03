# Level 1 Module DYN
# Simulates the dynamic behavior of the spacecraft
# Inputs: DYN initial conditions and its own outputs (for propagation)

from .DYN_TIME import DYN_TIME
from .DYN_ATT  import DYN_ATT
from .DYN_TRA  import DYN_TRA
from .DYN_SUN  import DYN_SUN

# Module output dictionary
DYN_out = {
    "DYN_TIME" : DYN_TIME.DYN_TIME_out,
    "DYN_ATT"  : DYN_ATT.DYN_ATT_out,
    "DYN_TRA"  : DYN_TRA.DYN_TRA_out,
    "DYN_SUN"  : DYN_SUN.DYN_SUN_out,
}

# Module main function
def run():
    DYN_TIME_out = DYN_TIME.run()
    DYN_ATT_out  = DYN_ATT.run()
    DYN_TRA_out  = DYN_TRA.run()
    DYN_SUN_out  = DYN_SUN.run()

    # Attribute outputs to DYN output
    DYN_out["DYN_TIME"] = DYN_TIME_out
    DYN_out["DYN_ATT"]  = DYN_ATT_out
    DYN_out["DYN_TRA"]  = DYN_TRA_out
    DYN_out["DYN_SUN"]  = DYN_SUN_out
    return dict(DYN_out)
