# DYN_TRA Module Test
# Use the main code solver to propagate spacecraft trajectory
# Compares the result using data from Spice
# Run from repository root with:
#   python3 -m DYN.DYN_TRA.Tests.test_DYN_TRA

import numpy as np
from DYN.DYN       import DYN
from PPC           import PPC
from PPC.PPC_plots import PPC_plots
from Utils         import spice
from Utils         import events
from Utils         import integrator

#########
# SETUP #
#########

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

##################
# MRO PARAMETERS #
##################

# Load NASA's MRO data with SPICE
# https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data
kernel_dir = "DYN/DYN_TRA/Tests/kernels/"
spice.load_kernel(kernel_dir + "mro_cruise.bsp")
spice.load_kernel(kernel_dir + "mro_sclkscet_00021_65536.tsc")

# Get Spacecraft initial condition
time_UTC_ini = spice.get_time("2005-08-15 T00:00:00")
MROpos_SSB_ini, MROvel_SSB_ini = spice.get_state("MRO", time_UTC_ini)

######################
# INITIALIZE MODULES #
######################

DYN_TIME_par = {
    "time_SIM_ini" : 0,
    "time_UTC_ini" : time_UTC_ini,
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
    DYN_obj.DYN_SUN,
    DYN_obj.DYN_EARTH,
    DYN_obj.DYN_MARS,
    DYN_obj.DYN_TRA,
    DYN_obj.DYN_GRV,
]

# Define MRO structure to store its states
MRO_TRA = {
    "MROpos_SSB" : MROpos_SSB_ini,
    "MROvel_SSB" : MROvel_SSB_ini,
}

# Initialize or load timeline
if DYN_log_load:
    DYN_timeline = PPC.load_timeline(DYN_log_path)
    timeline = {"DYN" : DYN_timeline}
else:
    timeline = {"DYN" : PPC.init_timeline(DYN_obj.snapshot(), n_steps),
                "MRO" : PPC.init_timeline(MRO_TRA, n_steps)}

########
# LOOP #
########

# Main loop
for step in range(n_steps):
    # Load external inputs
    inputs = events_table[sim_time]

    # Load DYN from file or compute everything
    if DYN_log_load:
        DYN_out = PPC.load_module(timeline["DYN"], step)
    else:
        # Update MRO state
        time_UTC = DYN_obj.snapshot()["DYN_TIME"]["time_UTC"]
        MROpos_SSB, MROvel_SSB = spice.get_state("MRO", time_UTC)
        # Update algebraic modules first
        DYN_obj.update_algebraic(sim_time, inputs)
        # Integrate all dynamic states together
        DYN_out = integrator.rk4_step(sim_time, sim_dt, DYN_obj, inputs)
        # Save results into the timeline
        output = {"DYN" : DYN_obj.snapshot(),
                  "MRO" : {"MROpos_SSB" : MROpos_SSB, "MROvel_SSB" : MROvel_SSB}}
        PPC.update_timeline(timeline, output, step)

    # Update time
    sim_time += sim_dt

############
# TEARDOWN #
############

# Save log (currently only for DYN module)
if DYN_log_save:
    PPC.store_timeline(timeline, DYN_log_path)

# Clear Kernels from memory
spice.clear_kernels()

###########
# RESULTS #
###########

# Plot spacecraft position vs MRO position
fig, ax_x = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0], ylabel="posx_SSB [km]", label="SC", title="MRO vs SC pos", subplot=(3,1,1))
fig, ax_y = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1], ylabel="posy_SSB [km]", label="SC", fig=fig, subplot=(3,1,2))
fig, ax_z = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2], xlabel="Time SIM [s]", ylabel="posz_SSB [km]", label="SC", fig=fig, subplot=(3,1,3))
PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["MRO"]["MROpos_SSB"][:,0], label="MRO", fig=fig, ax=ax_x)
PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["MRO"]["MROpos_SSB"][:,1], label="MRO", fig=fig, ax=ax_y)
PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["MRO"]["MROpos_SSB"][:,2], label="MRO", fig=fig, ax=ax_z)

# Plot spacecraft velocity vs MRO velocity
fig, ax_x = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCvel_SSB"][:,0], ylabel="velx_SSB [km/s]", label="SC", title="MRO vs SC vel", subplot=(3,1,1))
fig, ax_y = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCvel_SSB"][:,1], ylabel="vely_SSB [km/s]", label="SC", fig=fig, subplot=(3,1,2))
fig, ax_z = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCvel_SSB"][:,2], xlabel="Time SIM [s]", ylabel="velz_SSB [km/s]", label="SC", fig=fig, subplot=(3,1,3))
PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["MRO"]["MROvel_SSB"][:,0], label="MRO", fig=fig, ax=ax_x)
PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["MRO"]["MROvel_SSB"][:,1], label="MRO", fig=fig, ax=ax_y)
PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["MRO"]["MROvel_SSB"][:,2], label="MRO", fig=fig, ax=ax_z)

# Plot error on position comparing SC and MRO
SCpos_SSB  = timeline["DYN"]["DYN_TRA"]["SCpos_SSB"]
MROpos_SSB = timeline["MRO"]["MROpos_SSB"]
pos_error_x = SCpos_SSB[:,0] - MROpos_SSB[:,0]
pos_error_y = SCpos_SSB[:,1] - MROpos_SSB[:,1]
pos_error_z = SCpos_SSB[:,2] - MROpos_SSB[:,2]
pos_error = np.sqrt(pos_error_x**2 + pos_error_y**2 + pos_error_z**2)
fig, ax_x = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], pos_error_x, ylabel="pos_error_x [km]", label="SC", title="Error MRO vs SC pos", subplot=(4,1,1))
fig, ax_y = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], pos_error_y, ylabel="pos_error_y [km]", label="SC", fig=fig, subplot=(4,1,2))
fig, ax_z = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], pos_error_z, ylabel="pos_error_z [km]", label="SC", fig=fig, subplot=(4,1,3))
fig, ax_t = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], pos_error, xlabel="Time SIM [s]", ylabel="pos_error [km]", label="SC", fig=fig, subplot=(4,1,4))

# Plot error on velocity comparing SC and MRO
SCvel_SSB  = timeline["DYN"]["DYN_TRA"]["SCvel_SSB"]
MROvel_SSB = timeline["MRO"]["MROvel_SSB"]
vel_error_x = SCvel_SSB[:,0] - MROvel_SSB[:,0]
vel_error_y = SCvel_SSB[:,1] - MROvel_SSB[:,1]
vel_error_z = SCvel_SSB[:,2] - MROvel_SSB[:,2]
vel_error = np.sqrt(vel_error_x**2 + vel_error_y**2 + vel_error_z**2)
fig, ax_x = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], vel_error_x, ylabel="vel_error_x [km]", label="SC", title="Error MRO vs SC vel", subplot=(4,1,1))
fig, ax_y = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], vel_error_y, ylabel="vel_error_y [km]", label="SC", fig=fig, subplot=(4,1,2))
fig, ax_z = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], vel_error_z, ylabel="vel_error_z [km]", label="SC", fig=fig, subplot=(4,1,3))
fig, ax_t = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], vel_error, xlabel="Time SIM [s]", ylabel="vel_error [km]", label="SC", fig=fig, subplot=(4,1,4))

PPC.show_plot()