# DYN_TRA Module Test
# Use the main code solver to propagate spacecraft trajectory
# Compares the result using data from Spice
# Run from repository root with:
#   python3 -m DYN.DYN_TRA.Tests.test_DYN_TRA

import numpy as np
from DYN.DYN import DYN
from PPC     import PPC
from Utils   import spice
from Utils   import events
from Utils   import integrator

##################
# MRO PARAMETERS #
##################

# Load NASA's MRO data with SPICE
# https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data
kernel_dir = "Utils/kernels/"
spice.load_kernel(kernel_dir + "mro_cruise.bsp")
spice.load_kernel(kernel_dir + "mro_sclkscet_00021_65536.tsc")

# Get Spacecraft initial condition
# MRO trajectory will be compared during the cruise phase, in a period of time without any trajectory correction maneuvers
time_TDB_ini = spice.get_time("2005-10-01 T00:00:00")
time_TDB_end = spice.get_time("2005-10-31 T00:00:00")

MROpos_SSB_ini, MROvel_SSB_ini = spice.get_state("MRO", time_TDB_ini)

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
    "SCpos_SSB_ini" : MROpos_SSB_ini, # [km]   Initial position relative to SSB frame
    "SCvel_SSB_ini" : MROvel_SSB_ini, # [km/s] Initial velocity relative to SSB frame
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

# Define MRO structure to store its states
MRO_TRA = {
    "MROpos_SSB" : MROpos_SSB_ini,
    "MROvel_SSB" : MROvel_SSB_ini,
}

# Initialize timeline
timeline = {"DYN" : PPC.init_timeline(DYN_obj.snapshot(), n_steps),
            "MRO" : PPC.init_timeline(MRO_TRA, n_steps)}

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

    # Update MRO state
    time_TDB = DYN_obj.snapshot()["DYN_TIME"]["time_TDB"]
    MROpos_SSB, MROvel_SSB = spice.get_state("MRO", time_TDB)

    # Save results into the timeline
    output = {"DYN" : DYN_obj.snapshot(),
                "MRO" : {"MROpos_SSB" : MROpos_SSB, "MROvel_SSB" : MROvel_SSB}}
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

# Plot spacecraft position vs MRO position
fig, ax_x = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0], ylabel="posx_SSB [km]", label="SC", title="pos SC vs MRO", subplot=(3,1,1))
fig, ax_y = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1], ylabel="posy_SSB [km]", label="SC", fig=fig, subplot=(3,1,2))
fig, ax_z = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2], xlabel="Time SIM [days]", ylabel="posz_SSB [km]", label="SC", fig=fig, subplot=(3,1,3))
PPC.plot(time_days, timeline["MRO"]["MROpos_SSB"][:,0], label="MRO", fig=fig, ax=ax_x)
PPC.plot(time_days, timeline["MRO"]["MROpos_SSB"][:,1], label="MRO", fig=fig, ax=ax_y)
PPC.plot(time_days, timeline["MRO"]["MROpos_SSB"][:,2], label="MRO", fig=fig, ax=ax_z)

# Plot spacecraft velocity vs MRO velocity
fig, ax_x = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCvel_SSB"][:,0], ylabel="velx_SSB [km/s]", label="SC", title="vel SC vs MRO", subplot=(3,1,1))
fig, ax_y = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCvel_SSB"][:,1], ylabel="vely_SSB [km/s]", label="SC", fig=fig, subplot=(3,1,2))
fig, ax_z = PPC.plot(time_days, timeline["DYN"]["DYN_TRA"]["SCvel_SSB"][:,2], xlabel="Time SIM [days]", ylabel="velz_SSB [km/s]", label="SC", fig=fig, subplot=(3,1,3))
PPC.plot(time_days, timeline["MRO"]["MROvel_SSB"][:,0], label="MRO", fig=fig, ax=ax_x)
PPC.plot(time_days, timeline["MRO"]["MROvel_SSB"][:,1], label="MRO", fig=fig, ax=ax_y)
PPC.plot(time_days, timeline["MRO"]["MROvel_SSB"][:,2], label="MRO", fig=fig, ax=ax_z)

# Plot error on position comparing SC and MRO
SCpos_SSB  = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"]
MROpos_SSB = timeline["MRO"]["MROpos_SSB"]
pos_error_x = np.sqrt((SCpos_SSB[:,0] - MROpos_SSB[:,0])**2)
pos_error_y = np.sqrt((SCpos_SSB[:,1] - MROpos_SSB[:,1])**2)
pos_error_z = np.sqrt((SCpos_SSB[:,2] - MROpos_SSB[:,2])**2)
pos_error = np.sqrt(pos_error_x**2 + pos_error_y**2 + pos_error_z**2)
fig, ax_r = PPC.plot(time_days, pos_error_x, ylabel="pos error [km]", label="|$r_{SC,x} - r_{MRO,x}$|", title="Error pos SC vs MRO", subplot=(2,1,1), color=PPC.colors['blue'])
PPC.plot(time_days, pos_error_y, label="|$r_{SC,y} - r_{MRO,y}$|", fig=fig, ax=ax_r, color=PPC.colors['red'])
PPC.plot(time_days, pos_error_z, label="|$r_{SC,z} - r_{MRO,z}$|", fig=fig, ax=ax_r, color=PPC.colors['green'])
PPC.plot(time_days, pos_error, label="|$r_{SC} - r_{MRO}$|", fig=fig, ax=ax_r, color=PPC.colors['magenta'])

# Plot error on velocity comparing SC and MRO
SCvel_SSB  = timeline["DYN"]["DYN_TRA"]["SCvel_SSB"]
MROvel_SSB = timeline["MRO"]["MROvel_SSB"]
vel_error_x = np.sqrt((SCvel_SSB[:,0] - MROvel_SSB[:,0])**2)
vel_error_y = np.sqrt((SCvel_SSB[:,1] - MROvel_SSB[:,1])**2)
vel_error_z = np.sqrt((SCvel_SSB[:,2] - MROvel_SSB[:,2])**2)
vel_error = np.sqrt(vel_error_x**2 + vel_error_y**2 + vel_error_z**2)
fig, ax_v = PPC.plot(time_days, vel_error_x, ylabel="vel error [km/s]", label="|$v_{SC,x} - v_{MRO,x}$|", title="Error vel SC vs MRO", subplot=(2,1,2), color=PPC.colors['blue'], fig=fig)
PPC.plot(time_days, vel_error_y, label="|$v_{SC,y} - v_{MRO,y}$|", subplot=(4,1,2), fig=fig, ax=ax_v, color=PPC.colors['red'])
PPC.plot(time_days, vel_error_z, label="|$v_{SC,z} - v_{MRO,z}$|", subplot=(4,1,3), fig=fig, ax=ax_v, color=PPC.colors['green'])
PPC.plot(time_days, vel_error, xlabel="Time SIM [days]", ylabel="vel error [km/s]", label="|$v_{SC} - v_{MRO}$|", fig=fig, ax=ax_v, color=PPC.colors['magenta'])

# Trajectoy in the orbital plane
fig, ax = PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1], xlabel="x SSB [km]", ylabel="y SSB [km]", title="SC vs MRO trajectory", aspect='equal', color=PPC.colors["magenta"], style=':', zorder=4)
PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][-1,0], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][-1,1], label='SC', style='x', fig=fig, ax=ax, color=PPC.colors["magenta"], zorder=4)
PPC.plot(timeline["MRO"]["MROpos_SSB"][:,0],  timeline["MRO"]["MROpos_SSB"][:,1], fig=fig, ax=ax, color=PPC.colors["grey"])
PPC.plot(timeline["MRO"]["MROpos_SSB"][-1,0], timeline["MRO"]["MROpos_SSB"][-1,1], label="MRO", style='X', fig=fig, ax=ax, color=PPC.colors["grey"])
PPC.plot(timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"][:,0],  timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"][:,1],  fig=fig, ax=ax, color=PPC.colors["blue"])
PPC.plot(timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"][-1,0], timeline["DYN"]["DYN_EPH"]["EARTHpos_SSB"][-1,1], label="Earth", style='o', fig=fig, ax=ax, color=PPC.colors["blue"])
PPC.plot(timeline["DYN"]["DYN_EPH"]["MARSpos_SSB"][:,0],   timeline["DYN"]["DYN_EPH"]["MARSpos_SSB"][:,1],   fig=fig, ax=ax, color=PPC.colors["red"])
PPC.plot(timeline["DYN"]["DYN_EPH"]["MARSpos_SSB"][-1,0],  timeline["DYN"]["DYN_EPH"]["MARSpos_SSB"][-1,1],  label="Mars", style='o', fig=fig, ax=ax, color=PPC.colors["red"])
PPC.plot(timeline["DYN"]["DYN_EPH"]["SUNpos_SSB"][0,0],    timeline["DYN"]["DYN_EPH"]["SUNpos_SSB"][0,1],    label="Sun", style='o', fig=fig, ax=ax, color=PPC.colors["orange"])

PPC.show_plot()
