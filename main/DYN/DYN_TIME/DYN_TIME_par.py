from Utils import spice

# Pre computation of parameters
time_SIM_ini = 0 # [s] Initial simulation time
time_ET_ini = spice.get_time("2005-08-15 T00:00:00") # [s] ET initial time relative to J2000

# Parameters for Module DYN_TIME
DYN_TIME_par = {
    "time_SIM_ini" : time_SIM_ini,
    "time_ET_ini"  : time_ET_ini,
}
