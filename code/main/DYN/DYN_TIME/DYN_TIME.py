# Level 2 Module DYN_TIME
# Simulates the propagation of time for the simulation

from .DYN_TIME_par import DYN_TIME_par

# Module output dictionary
DYN_TIME_out = {
    "time" : DYN_TIME_par["time_ini"]
}

# Module main function
def run():
    DYN_TIME_out["time"] += DYN_TIME_par["dt"]
    return dict(DYN_TIME_out)
