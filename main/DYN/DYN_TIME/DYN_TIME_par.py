from Utils import spice

# Pre computation of parameters
dt           = 60 # [s] Simulation time step
time_SIM_ini = 0 # [s] Initial simulation time
# time_UTC_ini = spice.get_time("2025-10-14 17:00:00") # [s] UTC initial time relative to J2000
time_UTC_ini = spice.get_time("2005-08-15 T00:00:00") # [s] UTC initial time relative to J2000

# Parameters for Module DYN_TIME
DYN_TIME_par = {
    "dt"           : dt,
    "time_SIM_ini" : time_SIM_ini,
    "time_UTC_ini" : time_UTC_ini,
}
