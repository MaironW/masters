# Level 1 Module DYN
# Simulates the dynamic behavior of the spacecraft
# Inputs: DYN initial conditions and its own outputs (for propagation)

from .DYN_TIME import DYN_TIME

# Module output dictionary
DYN_out = {
    "DYN_TIME" : DYN_TIME.DYN_TIME_out
}

# Module main function
def run():
    DYN_TIME_out = DYN_TIME.run()

    # Attribute outputs to DYN output
    DYN_out["DYN_TIME"] = DYN_TIME_out
    return dict(DYN_out)
