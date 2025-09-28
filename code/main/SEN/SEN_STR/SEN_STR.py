# Level 2 Module SEN_STR
# Simulates a Star Tracker, giving spacecraft attitude and relative positions of selected selestial bodies
# Inputs: The spacecraft real attitude and position of celestial bodies

from .SEN_STR_par import SEN_STR_par

# Module output dictionary
SEN_STR_out = {
    "time" : SEN_STR_par["time_ini"]
}

# Module main function
def run(SEN_in):
    SEN_STR_out["time"] += SEN_STR_par["dt"]
    return dict(SEN_STR_out)
