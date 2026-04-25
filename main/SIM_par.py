# Simulation parameters

SIM_par = {
    # Time parameters
    "dt"           : 600,        # [s] 10 min
    "time_start"   : 0,          # [s]
    "time_end"     : 3600*24*10, # [s] 10 days

    # Log parameters
    "DYN_log_save" : False,
    "DYN_log_load" : False,
    "DYN_log_path" : "Logs/DYN",

    # List of default plots
    "PPC_plot_list" : [
        # "DYN_TIME",
        # "DYN_EPH",
        # "DYN_TRA",
        # "DYN_ATT",
        # "DYN_STR",
        # "DYN_PSR",
        # "SEN_TIME",
        # "SEN_STR",
        # "SEN_PSR",
        # "SEN_CMB",
        # "NAV_EPH",
        # "NAV_STR",
        # "NAV_CEL",
        # "NAV_PSR",
        # "NAV_CMB",
        # "NAV_EKF",
        # "NAV_UKF",
    ],
}
