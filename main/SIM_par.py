# Simulation parameters

SIM_par = {
    # Time parameters
    "dt"           : 600,       # [s] 10 min
    "time_start"   : 0,         # [s]
    "time_end"     : 3600*24*5, # [s] 5 days

    # Log parameters
    "DYN_log_save" : False,
    "DYN_log_load" : False,
    "DYN_log_path" : "Logs/DYN",

    # List of default plots
    "PPC_plot_list" : [
        # "DYN_TIME",
        # "DYN_SUN",
        # "DYN_EARTH",
        # "DYN_MARS",
        "DYN_TRA",
        # "DYN_GRV",
        # "DYN_ATT",
        # "DYN_STR",
        # "SEN_STR",
        # "NAV_CEL",
    ],
}
