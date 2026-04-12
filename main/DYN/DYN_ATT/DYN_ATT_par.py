import numpy as np

# Parameters for Module DYN_ATT

DYN_ATT_par = {
    "BOFq_SSB_ini" : np.array([1,0,0,0]), # Initial relative attitude from SSB to BOF reference frame
    "SCw_BOF_ini"  : np.array([0,0,0])    # [rad] Initial spacecraft angular speed on BOF frame
}
