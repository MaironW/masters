# DYN_EPH Module Test
# Based on test_DYN_TRA_MSL, just plots celestial bodies ephemerides from Spice
# Run from repository root with:
#   python3 -m DYN.DYN_EPH.Tests.test_DYN_EPH

import numpy as np
from DYN.DYN import DYN
from PPC     import PPC
from Utils   import spice
from Utils   import events
from Utils   import integrator


# Get Spacecraft initial condition
# MSL trajectory will be compared during the cruise phase, in a period of time without any trajectory correction maneuvers
time_TDB_ini = spice.get_time("2012-01-01 T00:00:00")
time_TDB_end = spice.get_time("2012-07-30 T00:00:00")

#########
# SETUP #
#########

SIM_par = {
    # Time parameters
    "dt"           : 600, # [s] 10 min
    "time_start"   : 0,   # [s]
    "time_end"     : int(time_TDB_end - time_TDB_ini), # [s]

    # Log parameters
    "DYN_log_save" : False,
    "DYN_log_load" : False,
    "DYN_log_path" : "Logs/DYN",

    # List of default plots
    "PPC_plot_list" : [
        "DYN_TIME",
        "DYN_EPH",
        "DYN_TRA",
    ],
}

# Load data or run new simulation
DYN_log_save = SIM_par["DYN_log_save"]
DYN_log_load = SIM_par["DYN_log_load"]
DYN_log_path = SIM_par["DYN_log_path"]

# Simulation parameters
sim_dt         = SIM_par["dt"]
sim_time_start = SIM_par["time_start"]
sim_time_end   = SIM_par["time_end"]
sim_time       = sim_time_start
step           = 0
n_steps        = int((sim_time_end - sim_time_start)/sim_dt) + 1

# Initialize inputs table
events_table = events.build_events_table(sim_time_start, sim_time_end, sim_dt)

######################
# INITIALIZE MODULES #
######################

DYN_TIME_par = {
    "time_SIM_ini" : 0,
    "time_TDB_ini"  : time_TDB_ini,
}

par_override = {
    "DYN_TIME" : DYN_TIME_par,
}

# Normally initialize DYN module, with override parameters
DYN_obj = DYN(par_override)
# Set only modules needed for the test
DYN_obj.modules = [
    DYN_obj.DYN_TIME,
    DYN_obj.DYN_EPH,
    DYN_obj.DYN_EPH,
    DYN_obj.DYN_EPH,
    DYN_obj.DYN_TRA,
]

# Initialize timeline
timeline = {"DYN" : PPC.init_timeline(DYN_obj.snapshot(), n_steps)}

########
# LOOP #
########

# Main loop
for step in range(1, n_steps):
    # Update time
    sim_time += sim_dt

    # Load external inputs
    inputs = events_table[sim_time]

    # Integrate all dynamic states together
    DYN_obj = integrator.rk4_step(sim_time, sim_dt, DYN_obj, inputs)
    # Update algebraic modules
    DYN_obj.update_algebraic(sim_time, None, inputs)

    # Update MSL state
    time_TDB = DYN_obj.snapshot()["DYN_TIME"]["time_TDB"]

    # Save results into the timeline
    output = {"DYN" : DYN_obj.snapshot()}
    PPC.update_timeline(timeline, output, step)

############
# TEARDOWN #
############

# Clear Kernels from memory
spice.clear_kernels()

###########
# RESULTS #
###########

EARTHpos_SSB  = timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"]
MOONpos_SSB   = timeline["DYN"]["DYN_EPH"]["MOONpos_SSB"]
MARSpos_SSB   = timeline["DYN"]["DYN_EPH"]["MARSpos_SSB"]
MARSvel_SSB   = timeline["DYN"]["DYN_EPH"]["MARSvel_SSB"]
DEIMOSpos_SSB = timeline["DYN"]["DYN_EPH"]["DEIMOSpos_SSB"]
PHOBOSpos_SSB = timeline["DYN"]["DYN_EPH"]["PHOBOSpos_SSB"]
SUNpos_SSB    = timeline["DYN"]["DYN_EPH"]["SUNpos_SSB"]

time_days = timeline["DYN"]["DYN_TIME"]["time_SIM"]/3600/24

# Solar System
fig, ax = PPC.plot([0], [0], [0], style='+', label="SSB", xlabel="X SSB [km]", ylabel="Y SSB [km]", zlabel="Z SSB [km]", title="Simulation Trajectories", aspect="equal",                      color=PPC.colors["black"])
bodies = [
    dict(name="SUN",    idx=0, color=PPC.colors["orange"]),
    dict(name="MOON",   idx=2, color=PPC.colors["grey"]),
    dict(name="EARTH",  idx=1, color=PPC.colors["blue"]),
    dict(name="DEIMOS", idx=4, color=PPC.colors["darkgrey"]),
    dict(name="PHOBOS", idx=5, color=PPC.colors["lightgrey"]),
    dict(name="MARS",   idx=3, color=PPC.colors["red"]),
]
for body in bodies:
    name = body["name"]
    PPC.plot(timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"][:,0], timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"][:,1], timeline["DYN"]["DYN_EPH"][f"{name}pos_SSB"][:,2], label=name, fig=fig, ax=ax, color=body["color"])

PPC.show_plot()

# Compute extra data, just to plot on text
MOONpos_ECI   = MOONpos_SSB - EARTHpos_SSB
DEIMOSpos_MCI = DEIMOSpos_SSB - MARSpos_SSB
PHOBOSpos_MCI = PHOBOSpos_SSB - MARSpos_SSB

# Save data
save_flg = True
if save_flg:
    # Downsample data
    n_samples = 200
    time_days_ds     = PPC.downsample(time_days,     n_samples)
    EARTHpos_SSB_ds  = PPC.downsample(EARTHpos_SSB,  n_samples)
    MARSpos_SSB_ds   = PPC.downsample(MARSpos_SSB,   n_samples)
    SUNpos_SSB_ds    = PPC.downsample(SUNpos_SSB,    n_samples)
    MOONpos_SSB_ds   = PPC.downsample(MOONpos_SSB,   n_samples)
    DEIMOSpos_SSB_ds = PPC.downsample(DEIMOSpos_SSB, n_samples)
    PHOBOSpos_SSB_ds = PPC.downsample(PHOBOSpos_SSB, n_samples)
    MOONpos_ECI_ds   = PPC.downsample(MOONpos_ECI,   n_samples)
    DEIMOSpos_MCI_ds = PPC.downsample(DEIMOSpos_MCI, n_samples)
    PHOBOSpos_MCI_ds = PPC.downsample(PHOBOSpos_MCI, n_samples)

    # Write data into file
    np.savetxt("EARTHpos_SSB.csv",  EARTHpos_SSB_ds,  delimiter=',')
    np.savetxt("MARSpos_SSB.csv",   MARSpos_SSB_ds,   delimiter=',')
    np.savetxt("SUNpos_SSB.csv",    SUNpos_SSB_ds,    delimiter=',')
    np.savetxt("MOONpos_SSB.csv",   MOONpos_SSB_ds,   delimiter=',')
    np.savetxt("DEIMOSpos_SSB.csv", DEIMOSpos_SSB_ds, delimiter=',')
    np.savetxt("PHOBOSpos_SSB.csv", PHOBOSpos_SSB_ds, delimiter=',')
