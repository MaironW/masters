# Level 2 Module DYN_ATT
# Simulates the propagation of the spacecraft attitude for the simulation

from .DYN_ATT_par import DYN_ATT_par

# Module output dictionary
def initialize():
    DYN_ATT_out = {
        "SSBq_BOF" : DYN_ATT_par["SSBq_BOF_ini"]
    }
    return DYN_ATT_out

# Module main function
def outputs(t, DYN_out, inputs):
    DYN_out["DYN_ATT"]["SSBq_BOF"] = inputs["DYN_ATT"]["SSBq_BOF"]
    return DYN_out
