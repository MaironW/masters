# DYN_TRA Module Test
# Use the main code solver to propagate spacecraft trajectory
# Checks if the trajectory is valid and going to Mars
# Run from repository root with:
#   python3 -m DYN.DYN_TRA.Tests.test_DYN_TRA_albuquerque

import numpy as np
from DYN.DYN       import DYN
from PPC           import PPC
from Utils         import spice
from Utils         import events
from Utils         import integrator
from Utils         import misc

# Spacecraft orbit configuration (From Albuquerque et al. 2024)
time_TDB_ini = spice.get_time("2023-04-15 T00:00:00")
time_TDB_end = spice.get_time("2023-11-04 T00:00:00")

# The orbit initial conditions are defined in ecliptic J2000,
# They need to be converted to standard (equatorial) J2000 coordinates before propagation
SCpos_SSB_ini_ECL = np.array([-1.3765e8, -6.2494e7, 3.2994e3]) # [km]   Initial position relative to SSB frame
SCvel_SSB_ini_ECL = np.array([-15.2727, -26.2104, -0.3666])    # [km/s] Initial velocity relative to SSB frame

SCpos_SSB_ini = misc.ECLtoEQT(SCpos_SSB_ini_ECL) # [km]
SCvel_SSB_ini = misc.ECLtoEQT(SCvel_SSB_ini_ECL) # [km/s]

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

DYN_TRA_par = {
    "ref_elements"  : "rvi",
    "BODY_ini"      : "SUN",
    "sma_ini"       : None,
    "ecc_ini"       : None,
    "incl_ini"      : None,
    "raan_ini"      : None,
    "argp_ini"      : None,
    "tano_ini"      : None,
    "SCpos_SSB_ini" : SCpos_SSB_ini, # [km]   Initial position relative to SSB frame
    "SCvel_SSB_ini" : SCvel_SSB_ini, # [km/s] Initial velocity relative to SSB frame
}

par_override = {
    "DYN_TIME" : DYN_TIME_par,
    "DYN_TRA"  : DYN_TRA_par,
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

time_days = timeline["DYN"]["DYN_TIME"]["time_SIM"]/3600/24

# Plot spacecraft position
fig, ax_x = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0], ylabel="posx_SSB [km]", label="SC", title="pos SC", subplot=(3,1,1))
fig, ax_y = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1], ylabel="posy_SSB [km]", label="SC", fig=fig, subplot=(3,1,2))
fig, ax_z = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2], xlabel="Time SIM [days]", ylabel="posz_SSB [km]", label="SC", fig=fig, subplot=(3,1,3))

# Plot spacecraft velocity
fig, ax_x = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCvel_SSB"][:,0], ylabel="velx_SSB [km/s]", label="SC", title="vel SC", subplot=(3,1,1))
fig, ax_y = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCvel_SSB"][:,1], ylabel="vely_SSB [km/s]", label="SC", fig=fig, subplot=(3,1,2))
fig, ax_z = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCvel_SSB"][:,2], xlabel="Time SIM [days]", ylabel="velz_SSB [km/s]", label="SC", fig=fig, subplot=(3,1,3))

# Trajectoy in the orbital plane
fig, ax = PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1], xlabel="x SSB [km]", ylabel="y SSB [km]", title="SC trajectory", aspect='equal', color=PPC.colors["magenta"], style=':', zorder=4)
PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][-1,0], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][-1,1], label='SC', style='x', fig=fig, ax=ax, color=PPC.colors["magenta"], zorder=4)
PPC.plot(timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"][:,0],  timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"][:,1],  fig=fig, ax=ax, color=PPC.colors["blue"])
PPC.plot(timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"][-1,0], timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"][-1,1], label="Earth", style='o', fig=fig, ax=ax, color=PPC.colors["blue"])
PPC.plot(timeline["DYN"]["DYN_EPH"]["MARSpos_SSB"][:,0],   timeline["DYN"]["DYN_EPH"]["MARSpos_SSB"][:,1],   fig=fig, ax=ax, color=PPC.colors["red"])
PPC.plot(timeline["DYN"]["DYN_EPH"]["MARSpos_SSB"][-1,0],  timeline["DYN"]["DYN_EPH"]["MARSpos_SSB"][-1,1],  label="Mars", style='o', fig=fig, ax=ax, color=PPC.colors["red"])
PPC.plot(timeline["DYN"]["DYN_EPH"]["SUNpos_SSB"][0,0],    timeline["DYN"]["DYN_EPH"]["SUNpos_SSB"][0,1],    label="Sun", style='o', fig=fig, ax=ax, color=PPC.colors["orange"])

# 3D trajectories
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
# Plot Spacecraft
PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2], label="SC", fig=fig, ax=ax, color=PPC.colors["magenta"])

PPC.show_plot()

