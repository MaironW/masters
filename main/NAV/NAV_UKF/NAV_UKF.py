# Level 2 Module NAV_UKF
# Provides the estimated spacecraft position and velocity with respect to the SSB frame
# Inputs: Calculated position, velocity and attitude from CELENAV, XNAV and CMBR

# Module output dictionary
NAV_UKF_out = {
    "dummy" : None
}

# Module main function
def run(NAV_CELENAV_out, NAV_XNAV_out, NAV_CMBR_out):
    return dict(NAV_UKF_out)
