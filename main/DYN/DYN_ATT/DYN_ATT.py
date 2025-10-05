# Level 2 Module DYN_ATT
# Simulates the propagation of the spacecraft attitude for the simulation

from .DYN_ATT_par import DYN_ATT_par

# Module output dictionary
DYN_ATT_out = {
    "SSBq_BOF" : DYN_ATT_par["SSBq_BOF_ini"]
}

# Module main function
def run(inputs):
    DYN_ATT_out["SSBq_BOF"] = inputs["SSBq_BOF"]
    return dict(DYN_ATT_out)
